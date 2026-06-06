"""
vm.py
=====

Execution layer for Causal Defect Algebra programs.

The DefectVM is intentionally small and backend-agnostic: it wraps a
DefectState or CausalDefectEngine, registers operations, and executes a list of
instructions.  This turns the package from analyzers/backends into a programmable
computational medium over defect states.

Instruction formats
-------------------
All of the following are accepted:

    {"op": "observe", "args": {"name": "energy"}, "save_as": "E"}
    ("observe", {"name": "energy"})
    "summary"

Programs are ordinary lists of instructions.  Operation outputs are stored in
``vm.context`` when ``save_as`` is provided.  Later args may reference prior
outputs with strings of the form ``"$name"``.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple, Union

import numpy as np

from .core import CausalDefectEngine, DefectState, Defect

Instruction = Union[str, Tuple[str, Mapping[str, Any]], Mapping[str, Any]]
Operation = Callable[..., Any]


@dataclass
class VMTraceEntry:
    """A single executed instruction and its result metadata."""

    op: str
    args: Dict[str, Any] = field(default_factory=dict)
    save_as: Optional[str] = None
    result_type: str = "None"
    error: Optional[str] = None


class DefectVM:
    """Programmable runtime for a fitted defect state or engine.

    Parameters
    ----------
    state_or_engine:
        Either a :class:`DefectState` or a fitted :class:`CausalDefectEngine`.
        If an engine is supplied, backend methods such as ``compose`` and
        ``holonomy`` are automatically exposed when available.
    strict:
        If True, unknown operations or operation failures raise immediately. If
        False, failures are written into the trace and execution continues.
    copy_state:
        Shallow-copy the state reference.  The state payload itself is not deep
        copied because it can contain large arrays.
    """

    def __init__(self, state_or_engine: Union[DefectState, CausalDefectEngine], *, strict: bool = True, copy_state: bool = False):
        self.engine: Optional[CausalDefectEngine]
        if isinstance(state_or_engine, CausalDefectEngine):
            self.engine = state_or_engine
            self.state = state_or_engine.state
        elif isinstance(state_or_engine, DefectState):
            self.engine = None
            self.state = state_or_engine
        else:
            raise TypeError("DefectVM expects a DefectState or CausalDefectEngine")

        if copy_state:
            self.state = DefectState(
                events=self.state.events,
                generators=self.state.generators,
                connection=self.state.connection,
                fusion_rules=self.state.fusion_rules,
                projectors=self.state.projectors,
                omega=self.state.omega,
                defects=list(self.state.defects),
                observables=dict(self.state.observables),
                metadata=dict(self.state.metadata),
                valid=self.state.valid,
            )

        self.strict = bool(strict)
        self.ops: Dict[str, Operation] = {}
        self.context: Dict[str, Any] = {"state": self.state}
        self.trace: List[VMTraceEntry] = []
        self.register_builtin_ops()

    # ------------------------------------------------------------------
    # Registration and program execution
    # ------------------------------------------------------------------

    def register(self, name: str, fn: Operation, *, overwrite: bool = True) -> "DefectVM":
        key = str(name).strip()
        if not key:
            raise ValueError("operation name cannot be empty")
        if not overwrite and key in self.ops:
            raise KeyError(f"operation already registered: {key}")
        self.ops[key] = fn
        return self

    def register_many(self, ops: Mapping[str, Operation], *, overwrite: bool = True) -> "DefectVM":
        for name, fn in ops.items():
            self.register(name, fn, overwrite=overwrite)
        return self

    def parse_instruction(self, instr: Instruction) -> Tuple[str, Dict[str, Any], Optional[str]]:
        if isinstance(instr, str):
            return instr, {}, None
        if isinstance(instr, tuple):
            if len(instr) == 2:
                op, args = instr
                return str(op), dict(args), None
            if len(instr) == 3:
                op, args, save_as = instr
                return str(op), dict(args), None if save_as is None else str(save_as)
            raise ValueError("tuple instructions must be (op,args) or (op,args,save_as)")
        if isinstance(instr, Mapping):
            if "op" not in instr:
                raise ValueError("dict instruction requires an 'op' key")
            op = str(instr["op"])
            args = dict(instr.get("args", {}))
            save_as = instr.get("save_as", instr.get("as", None))
            return op, args, None if save_as is None else str(save_as)
        raise TypeError(f"unsupported instruction type: {type(instr)!r}")

    def _resolve_value(self, value: Any) -> Any:
        if isinstance(value, str) and value.startswith("$"):
            key = value[1:]
            if key not in self.context:
                raise KeyError(f"context variable not found: {value}")
            return self.context[key]
        if isinstance(value, list):
            return [self._resolve_value(v) for v in value]
        if isinstance(value, tuple):
            return tuple(self._resolve_value(v) for v in value)
        if isinstance(value, dict):
            return {k: self._resolve_value(v) for k, v in value.items()}
        return value

    def _resolve_args(self, args: Mapping[str, Any]) -> Dict[str, Any]:
        return {k: self._resolve_value(v) for k, v in args.items()}

    def step(self, instr: Instruction) -> Any:
        op, raw_args, save_as = self.parse_instruction(instr)
        args = self._resolve_args(raw_args)
        entry = VMTraceEntry(op=op, args=dict(raw_args), save_as=save_as)
        try:
            if op not in self.ops:
                raise KeyError(f"unknown VM operation: {op}")
            result = self.ops[op](**args)
            entry.result_type = type(result).__name__
            if save_as is not None:
                self.context[save_as] = result
            self.context["last"] = result
            self.trace.append(entry)
            return result
        except Exception as exc:
            entry.error = str(exc)
            self.trace.append(entry)
            if self.strict:
                raise
            return None

    def run(self, program: Sequence[Instruction], *, reset_trace: bool = True) -> Any:
        if reset_trace:
            self.trace.clear()
        out = None
        for instr in program:
            out = self.step(instr)
        return out

    # ------------------------------------------------------------------
    # Built-in operations
    # ------------------------------------------------------------------

    def register_builtin_ops(self) -> None:
        self.register_many(
            {
                "summary": self.op_summary,
                "state": self.op_state,
                "observe": self.op_observe,
                "set_observable": self.op_set_observable,
                "defects": self.op_defects,
                "defect_table": self.op_defect_table,
                "compose": self.op_compose,
                "invert": self.op_invert,
                "power": self.op_power,
                "project": self.op_project,
                "holonomy": self.op_holonomy,
                "transport": self.op_holonomy,
                "omega": self.op_omega,
                "sector_norms": self.op_sector_norms,
                "select_defects": self.op_select_defects,
                "map": self.op_map,
            }
        )

    def op_state(self) -> DefectState:
        return self.state

    def op_summary(self) -> Dict[str, Any]:
        out = self.state.summary()
        out["context_keys"] = sorted(self.context.keys())
        out["ops"] = sorted(self.ops.keys())
        return out

    def op_observe(self, name: Optional[str] = None, default: Any = None) -> Any:
        if name is None:
            return dict(self.state.observables)
        return self.state.observables.get(name, default)

    def op_set_observable(self, name: str, value: Any) -> Any:
        self.state.observables[str(name)] = value
        return value

    def op_defects(self, kind: Optional[str] = None, min_strength: Optional[float] = None, min_abs_charge: Optional[float] = None) -> List[Defect]:
        out = list(self.state.defects)
        if kind is not None:
            out = [d for d in out if d.kind == kind]
        if min_strength is not None:
            out = [d for d in out if d.strength >= float(min_strength)]
        if min_abs_charge is not None:
            out = [d for d in out if abs(d.charge) >= float(min_abs_charge)]
        return out

    def op_defect_table(self, **filters: Any) -> List[Dict[str, Any]]:
        return [d.as_dict() for d in self.op_defects(**filters)]

    def op_compose(self, a: Any, b: Any) -> Any:
        if self.state.fusion_rules is not None:
            return self.state.fusion_rules.compose(a, b)
        if self.engine is not None:
            return self.engine.compose(a, b)
        raise RuntimeError("no composition rule installed")

    def op_invert(self, a: Any) -> Any:
        if self.state.fusion_rules is None:
            raise RuntimeError("no composition rule installed")
        return self.state.fusion_rules.invert(a)

    def op_power(self, a: Any, k: int) -> Any:
        if self.state.fusion_rules is None:
            raise RuntimeError("no composition rule installed")
        return self.state.fusion_rules.power(a, int(k))

    def op_project(self, sector: Any, x: Any = None) -> Any:
        if self.state.projectors is None:
            raise RuntimeError("no projectors installed")
        if sector not in self.state.projectors:
            raise KeyError(f"unknown projector sector: {sector}")
        projector = self.state.projectors[sector]
        if x is None:
            x = self.state.generators
        return projector(x)

    def op_holonomy(self, path: Sequence[int], **kwargs: Any) -> Any:
        if self.engine is not None:
            return self.engine.holonomy(path, **kwargs) if kwargs else self.engine.holonomy(path)
        connection = self.state.connection
        if connection is None:
            raise RuntimeError("no connection installed")
        # Generic connection forms: dict[(i,j)] -> matrix/scalar, or square array.
        if isinstance(connection, Mapping):
            acc = None
            for a, b in zip(path[:-1], path[1:]):
                u = connection.get((int(a), int(b)))
                if u is None:
                    raise KeyError(f"no connection edge {(a,b)}")
                acc = u if acc is None else _compose_transport(acc, u)
            return acc
        arr = np.asarray(connection)
        vals = [arr[int(a), int(b)] for a, b in zip(path[:-1], path[1:])]
        if not vals:
            return 0.0
        return np.sum(vals, axis=0)

    def op_omega(self) -> Any:
        return self.state.omega

    def op_sector_norms(self) -> Dict[str, float]:
        """Return basic norms for generators/projectors/omega for quick diagnostics."""
        out: Dict[str, float] = {}
        if self.state.generators is not None:
            try:
                out["generators"] = float(np.linalg.norm(np.asarray(self.state.generators)))
            except Exception:
                pass
        if self.state.omega is not None:
            try:
                out["omega"] = float(np.linalg.norm(np.asarray(self.state.omega)))
            except Exception:
                pass
        if self.state.defects:
            out["defect_charge_l1"] = float(sum(abs(d.charge) for d in self.state.defects))
            out["defect_strength_l1"] = float(sum(abs(d.strength) for d in self.state.defects))
        return out

    def op_select_defects(self, key: str = "strength", top_k: Optional[int] = None, reverse: bool = True) -> List[Defect]:
        def val(d: Defect) -> float:
            if key == "strength":
                return float(d.strength)
            if key == "charge":
                return float(d.charge)
            if key == "abs_charge":
                return abs(float(d.charge))
            return float(d.metadata.get(key, 0.0))

        defects = sorted(self.state.defects, key=val, reverse=reverse)
        return defects if top_k is None else defects[: int(top_k)]

    def op_map(self, fn: Callable[[Any], Any], values: Sequence[Any]) -> List[Any]:
        return [fn(v) for v in values]


def _compose_transport(a: Any, b: Any) -> Any:
    """Generic transport composition helper."""
    aa = np.asarray(a)
    bb = np.asarray(b)
    if aa.ndim == 2 and bb.ndim == 2 and aa.shape[1] == bb.shape[0]:
        return bb @ aa
    return aa + bb
