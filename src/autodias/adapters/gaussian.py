from __future__ import annotations

from pathlib import Path

from autodias.adapters.base import EngineAdapter


class GaussianAdapter(EngineAdapter):
    name = "gaussian"

    def render_input(
        self,
        template: str,
        coordinates: str,
        charge: int,
        multiplicity: int,
        job_name: str,
    ) -> str:
        rendered = template.replace("$coordinates", coordinates.strip())
        rendered = rendered.replace("$charge", str(charge))
        rendered = rendered.replace("$multiplicity", str(multiplicity))
        rendered = rendered.replace("$filename", job_name)
        return rendered

    def parse_energy(self, output_path: Path) -> float:
        energy = 0.0
        for line in reversed(output_path.read_text(errors="replace").splitlines()):
            if "SCF Done" in line:
                energy = float(line.split()[4])
                break
            if "Energy=" in line:
                energy = float(line.split()[3])
                break
        return energy
