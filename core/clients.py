from __future__ import annotations
import socket
import time
from typing import Any, Dict, Iterable, List, Optional
from urllib.parse import urlparse
import requests


class BaseNIMClient:

    def __init__(
        self, base_url: str, timeout_seconds: int = 120, max_retries: int = 3
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.timeout_seconds = timeout_seconds
        self.max_retries = max_retries
        self.health_check_ttl_seconds = 15.0
        self._last_health_check_at = 0.0
        self._last_health_check_ok = False

    def _http_get_ok(self, path: str) -> bool:
        try:
            response = requests.get(
                f"{self.base_url}{path}", timeout=min(3, self.timeout_seconds)
            )
            return response.status_code < 500
        except requests.RequestException:
            return False

    def _tcp_port_alive(self) -> bool:
        parsed = urlparse(self.base_url)
        host = parsed.hostname
        port = parsed.port
        if not host or not port:
            return False
        try:
            with socket.create_connection((host, port), timeout=2):
                return True
        except OSError:
            return False

    def health_check(self, force: bool = False) -> bool:
        now = time.monotonic()
        if (
            not force
            and now - self._last_health_check_at <= self.health_check_ttl_seconds
        ):
            return self._last_health_check_ok
        if not self._tcp_port_alive():
            self._last_health_check_at = now
            self._last_health_check_ok = False
            return False
        for health_path in ("/health", "/v1/health", "/ready", "/"):
            if self._http_get_ok(health_path):
                self._last_health_check_at = now
                self._last_health_check_ok = True
                return True
        self._last_health_check_at = now
        self._last_health_check_ok = False
        return False

    def wait_for_healthy(
        self, max_wait_seconds: int = 300, poll_interval: int = 5
    ) -> bool:
        deadline = time.monotonic() + max_wait_seconds
        attempt = 0
        while time.monotonic() < deadline:
            attempt += 1
            if self.health_check(force=True):
                return True
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            time.sleep(min(poll_interval, remaining))
        return False


class Evo2Client(BaseNIMClient):

    def __init__(
        self,
        base_url: str = "http://localhost:8000",
        timeout_seconds: int = 120,
        max_retries: int = 3,
    ) -> None:
        super().__init__(
            base_url=base_url, timeout_seconds=timeout_seconds, max_retries=max_retries
        )

    @staticmethod
    def _extract_score(payload: Dict[str, Any]) -> float:
        for key in ("likelihood", "log_likelihood", "score", "probability"):
            value = payload.get(key)
            if isinstance(value, (int, float)):
                return float(value)
        data = payload.get("data")
        if isinstance(data, dict):
            for key in ("likelihood", "log_likelihood", "score", "probability"):
                value = data.get(key)
                if isinstance(value, (int, float)):
                    return float(value)
        raise ValueError("Unable to parse Evo2 likelihood score from response.")

    def get_likelihood(self, sequence: str) -> float:
        if not self.health_check():
            raise RuntimeError(f"Evo2 service not healthy at {self.base_url}")
        payload = {"sequence": sequence}
        last_error = None
        for attempt in range(self.max_retries):
            for endpoint in ("/likelihood", "/v1/likelihood", "/score", "/v1/score"):
                try:
                    response = requests.post(
                        f"{self.base_url}{endpoint}",
                        json=payload,
                        timeout=self.timeout_seconds,
                    )
                    response.raise_for_status()
                    return self._extract_score(response.json())
                except requests.RequestException as exc:
                    last_error = f"{endpoint}: {exc}"
            if attempt < self.max_retries - 1:
                time.sleep(2**attempt)
        raise RuntimeError(
            f"Evo2 failed after {self.max_retries} retries. Last error: {last_error}"
        )

    def get_delta(self, wt_sequence: str, mutant_sequence: str) -> float:
        wt_score = self.get_likelihood(wt_sequence)
        mt_score = self.get_likelihood(mutant_sequence)
        return mt_score - wt_score


def _flatten_numeric_values(data: Any, parent_key: str = "") -> Dict[str, float]:
    flattened: Dict[str, float] = {}
    if isinstance(data, dict):
        for key, value in data.items():
            new_key = f"{parent_key}.{key}" if parent_key else key
            flattened.update(_flatten_numeric_values(value, new_key))
    elif isinstance(data, list):
        for idx, item in enumerate(data):
            new_key = f"{parent_key}[{idx}]"
            flattened.update(_flatten_numeric_values(item, new_key))
    elif isinstance(data, (int, float)):
        flattened[parent_key] = float(data)
    return flattened


def parse_boltz2_metrics(response_json: Dict[str, Any]) -> Dict[str, Any]:
    values = _flatten_numeric_values(response_json)

    def first_match(candidates: Iterable[str]) -> Optional[float]:
        lower_map = {key.lower(): value for key, value in values.items()}
        for candidate in candidates:
            candidate_lower = candidate.lower()
            for key, value in lower_map.items():
                if candidate_lower in key:
                    return value
        return None

    return {
        "iptm": first_match(("iptm", "interface_tm", "interface_confidence")),
        "ptm": first_match(("ptm",)),
        "plddt": first_match(("plddt", "confidence", "confidence_score")),
        "complex_energy": first_match(
            ("complex_energy", "total_energy", "system_energy")
        ),
        "binding_affinity": first_match(("binding_affinity", "affinity")),
        "affinity_pred_value": first_match(("affinity_pred_value",)),
        "pair_energy": first_match(("pair_energy", "delta_g")),
        "stability_score": first_match(("stability", "folding_energy")),
        "clash_score": first_match(("clash", "clash_score", "steric_clash")),
        "artifact_dir": response_json.get("artifact_dir"),
    }


class Boltz2Client(BaseNIMClient):

    def __init__(
        self,
        base_url: str = "http://localhost:8001",
        timeout_seconds: int = 300,
        max_retries: int = 3,
    ) -> None:
        super().__init__(
            base_url=base_url, timeout_seconds=timeout_seconds, max_retries=max_retries
        )

    def _candidate_predict_endpoints(self) -> List[str]:
        if self.base_url.startswith("http://localhost:") or self.base_url.startswith(
            "http://127.0.0.1:"
        ):
            return ["/predict"]
        if str(
            __import__("os").getenv("AMR_BOLTZ2_FALLBACK_ENDPOINTS", "")
        ).strip().lower() in {"1", "true", "yes", "on"}:
            return ["/predict", "/v1/predict", "/structure", "/v1/structure"]
        return ["/predict"]

    def _candidate_prefetch_endpoints(self) -> List[str]:
        if self.base_url.startswith("http://localhost:") or self.base_url.startswith(
            "http://127.0.0.1:"
        ):
            return ["/prefetch"]
        if str(
            __import__("os").getenv("AMR_BOLTZ2_FALLBACK_ENDPOINTS", "")
        ).strip().lower() in {"1", "true", "yes", "on"}:
            return ["/prefetch", "/v1/prefetch", "/prepare", "/v1/prepare"]
        return ["/prefetch"]

    def _post_json(
        self, payload: Dict[str, Any], endpoints: List[str]
    ) -> Dict[str, Any]:
        last_error = None
        for attempt in range(self.max_retries):
            for endpoint in endpoints:
                try:
                    response = requests.post(
                        f"{self.base_url}{endpoint}",
                        json=payload,
                        timeout=self.timeout_seconds,
                    )
                    response.raise_for_status()
                    return response.json()
                except requests.RequestException as exc:
                    last_error = f"{endpoint}: {exc}"
            if attempt < self.max_retries - 1:
                time.sleep(2**attempt)
        raise RuntimeError(
            f"Boltz2 failed after {self.max_retries} retries. Last error: {last_error}"
        )

    def predict(
        self,
        scenario: str,
        protein_sequence: str,
        protein_pdb: Optional[str] = None,
        ligand_sdf: Optional[str] = None,
        dna_sequence: Optional[str] = None,
        partner_sequences: Optional[List[str]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        if not self.health_check():
            raise RuntimeError(f"Boltz2 service not healthy at {self.base_url}")
        payload = self._build_payload(
            scenario,
            protein_sequence,
            protein_pdb,
            ligand_sdf,
            dna_sequence,
            partner_sequences=partner_sequences,
        )
        if metadata:
            payload["_meta"] = metadata
        response_json = self._post_json(payload, self._candidate_predict_endpoints())
        return parse_boltz2_metrics(response_json)

    def prefetch(
        self,
        scenario: str,
        protein_sequence: str,
        protein_pdb: Optional[str] = None,
        ligand_sdf: Optional[str] = None,
        dna_sequence: Optional[str] = None,
        partner_sequences: Optional[List[str]] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        if not self.health_check():
            raise RuntimeError(f"Boltz2 service not healthy at {self.base_url}")
        payload = self._build_payload(
            scenario,
            protein_sequence,
            protein_pdb,
            ligand_sdf,
            dna_sequence,
            partner_sequences=partner_sequences,
        )
        if metadata:
            payload["_meta"] = metadata
        return self._post_json(payload, self._candidate_prefetch_endpoints())

    def _build_ligand_payload(
        self,
        protein_sequence: str,
        protein_pdb: Optional[str],
        ligand_sdf: Optional[str],
    ) -> Dict[str, Any]:
        return {
            "scenario": "protein_ligand",
            "protein": {"sequence": protein_sequence, "template_pdb": protein_pdb},
            "ligand": {"sdf_path": ligand_sdf},
        }

    def _build_dna_payload(
        self,
        protein_sequence: str,
        protein_pdb: Optional[str],
        dna_sequence: Optional[str],
        dimer: bool = False,
    ) -> Dict[str, Any]:
        if dimer:
            return {
                "scenario": "protein_dna",
                "chains": [
                    {"sequence": protein_sequence, "type": "protein"},
                    {"sequence": protein_sequence, "type": "protein"},
                ],
                "template_pdb": protein_pdb,
                "dna": {"sequence": dna_sequence},
            }
        return {
            "scenario": "protein_dna",
            "protein": {"sequence": protein_sequence, "template_pdb": protein_pdb},
            "dna": {"sequence": dna_sequence},
        }

    def _build_protein_payload(
        self, protein_sequence: str, protein_pdb: Optional[str]
    ) -> Dict[str, Any]:
        return {
            "scenario": "protein_only",
            "protein": {"sequence": protein_sequence, "template_pdb": protein_pdb},
        }

    def _build_complex_payload(
        self,
        protein_sequence: str,
        protein_pdb: Optional[str],
        partner_sequences: Optional[List[str]],
        dimer: bool = False,
    ) -> Dict[str, Any]:
        chains = [{"sequence": protein_sequence, "type": "protein"}]
        for partner_sequence in partner_sequences or []:
            if partner_sequence:
                chains.append({"sequence": partner_sequence, "type": "protein"})
        if dimer and len(chains) == 1:
            chains.append({"sequence": protein_sequence, "type": "protein"})
        if len(chains) < 2:
            return self._build_protein_payload(protein_sequence, protein_pdb)
        return {
            "scenario": "protein_protein",
            "chains": chains,
            "template_pdb": protein_pdb,
        }

    def _build_payload(
        self,
        scenario: str,
        protein_sequence: str,
        protein_pdb: Optional[str],
        ligand_sdf: Optional[str],
        dna_sequence: Optional[str],
        partner_sequences: Optional[List[str]] = None,
    ) -> Dict[str, Any]:
        normalized = (scenario or "").upper()
        if normalized == "PROTEIN_LIGAND":
            return self._build_ligand_payload(protein_sequence, protein_pdb, ligand_sdf)
        if normalized in ("PROTEIN_DNA", "DIMER_DNA"):
            return self._build_dna_payload(
                protein_sequence,
                protein_pdb,
                dna_sequence,
                dimer=normalized == "DIMER_DNA",
            )
        if normalized in ("PROTEIN_PROTEIN", "COMPLEX", "DIMER"):
            return self._build_complex_payload(
                protein_sequence,
                protein_pdb,
                partner_sequences,
                dimer=normalized == "DIMER",
            )
        if normalized in ("FOLDING_ONLY", "PROTEIN_ONLY"):
            return self._build_protein_payload(protein_sequence, protein_pdb)
        raise ValueError(f"Unsupported scenario: {scenario}")
