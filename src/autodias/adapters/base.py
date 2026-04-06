from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path


class EngineAdapter(ABC):
    name: str

    @abstractmethod
    def render_input(
        self,
        template: str,
        coordinates: str,
        charge: int,
        multiplicity: int,
        job_name: str,
    ) -> str:
        raise NotImplementedError

    @abstractmethod
    def parse_energy(self, output_path: Path) -> float:
        raise NotImplementedError
