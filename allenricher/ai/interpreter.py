"""Optional AI-assisted interpretation for AllEnricher result tables.

All backends implement the same interface. Statistical analysis remains in
AllEnricher; an interpreter only summarizes already-computed results."""

import os
import json
import logging
import math
import re
from abc import ABC, abstractmethod
from itertools import combinations
from typing import Dict, List, Optional, Any
import pandas as pd

logger = logging.getLogger(__name__)

AI_MODES = {"summary", "reviewer", "caption"}
_AI_SECTIONS = ("core_themes", "key_evidence", "limitations", "computational_checks")
_DEFAULT_EVIDENCE_LIMITS = {"ora": 15, "gsea": 10, "ssgsea": 10, "gsva": 10}


def _first_value(row: pd.Series, names: tuple[str, ...], default: Any = None) -> Any:
    """Return the first non-empty value from a row using known aliases."""
    for name in names:
        if name in row.index:
            value = row.get(name)
            if value is not None and not (isinstance(value, float) and math.isnan(value)):
                if str(value).strip().lower() not in {"", "nan", "none", "na"}:
                    return value
    return default


def _number(value: Any, default: Optional[float] = None) -> Optional[float]:
    try:
        number = float(value)
        return number if math.isfinite(number) else default
    except (TypeError, ValueError):
        return default


def _json_value(value: Any) -> Any:
    """Convert pandas/numpy scalar values into JSON-safe values."""
    if value is None:
        return None
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_json_value(item) for item in value]
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try:
            return _json_value(value.item())
        except (TypeError, ValueError):
            pass
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _split_genes(value: Any) -> List[str]:
    if value is None:
        return []
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "na"}:
        return []
    text = text.strip("[]()")
    return [gene.strip().strip("'\"") for gene in re.split(r"[;,/|\s]+", text) if gene.strip()]


def _term_id(row: pd.Series, fallback: Any) -> str:
    return str(_first_value(row, ("Term_ID", "ID", "term_id", "pathway", "Term"), fallback))


def _term_name(row: pd.Series, fallback: Any) -> str:
    return str(_first_value(row, ("Term_Name", "Description", "term_name", "Name", "pathway", "Term"), fallback))


def _database_prefix(database: str) -> str:
    known = (
        "PUBLIC_GO_CUSTOM", "WikiPathways", "AnimalTFDB", "DisGeNET",
        "hTFtarget", "Reactome", "TRRUST", "ChEA3", "KEGG", "CUSTOM", "DO", "GO",
    )
    raw = str(database)
    for prefix in known:
        if raw == prefix or raw.startswith(prefix + "_"):
            return prefix
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", raw).strip("._") or "Database"


def _is_tf_database(database: str, frame: pd.DataFrame) -> bool:
    token = str(database).lower()
    return any(name in token for name in ("trrust", "chea3", "animaltfdb", "htftarget", "tf")) or any(
        column in frame.columns for column in ("TF", "Target_Genes", "Target_Count", "Library", "Context")
    )


def _row_record(
    evidence_id: str,
    database: str,
    method: str,
    kind: str,
    position: int,
    row: pd.Series,
    values: Dict[str, Any],
) -> Dict[str, Any]:
    term_id = _term_id(row, position)
    term_name = _term_name(row, term_id)
    return {
        "evidence_id": evidence_id,
        "database": str(database),
        "method": method,
        "kind": kind,
        "row_position": position,
        "term_id": term_id,
        "term_name": term_name,
        "values": {key: _json_value(value) for key, value in values.items()},
        "raw": {str(key): _json_value(value) for key, value in row.to_dict().items()},
    }


def _rank_rows(frame: pd.DataFrame, q_names: tuple[str, ...], strength_names: tuple[str, ...]) -> List[tuple[int, pd.Series]]:
    ranked = []
    for position, (_, row) in enumerate(frame.iterrows()):
        q_value = _number(_first_value(row, q_names), 1.0)
        strength = abs(_number(_first_value(row, strength_names), 0.0) or 0.0)
        ranked.append((q_value if q_value is not None else 1.0, -strength, position, row))
    ranked.sort(key=lambda item: (item[0], item[1], item[2]))
    return [(item[2], item[3]) for item in ranked]


def _activity_sample_columns(frame: pd.DataFrame) -> List[str]:
    metadata = {
        "Term_ID", "Term_Name", "Description", "Hierarchy", "Database", "Gene_Count",
        "Background_Count", "Term_URL", "NES", "ES", "P_Value", "Adjusted_P_Value",
        "FDR", "Genes", "Leading_Edge", "leadingEdge", "Rich_Factor", "Gene_Ratio",
        "Background_Ratio", "Expected_Count", "pval", "padj", "size", "setSize",
        "pathway", "log2err",
    }
    return [
        str(column) for column in frame.columns
        if column not in metadata and pd.api.types.is_numeric_dtype(frame[column])
    ]


def _activity_summary(row: pd.Series, sample_columns: List[str], groups: Optional[Dict[str, List[str]]]) -> Dict[str, Any]:
    values = {sample: _number(row.get(sample)) for sample in sample_columns}
    group_means: Dict[str, float] = {}
    if groups:
        for group, samples in groups.items():
            observed = [values[sample] for sample in samples if sample in values and values[sample] is not None]
            if observed:
                group_means[group] = sum(observed) / len(observed)
    differences = {}
    for left, right in combinations(group_means, 2):
        differences[f"{left} vs {right}"] = group_means[left] - group_means[right]

    observed_values = [value for value in values.values() if value is not None]
    outliers: List[str] = []
    if len(observed_values) >= 4:
        series = pd.Series(observed_values)
        lower = series.quantile(0.25)
        upper = series.quantile(0.75)
        spread = upper - lower
        for sample, value in values.items():
            if value is not None and (value < lower - 1.5 * spread or value > upper + 1.5 * spread):
                outliers.append(sample)
    return {
        "sample_values": values,
        "group_means": group_means,
        "group_differences": differences,
        "outlier_samples": outliers,
    }


def build_structured_evidence(
    results: Dict[str, pd.DataFrame],
    method: str = "hypergeometric",
    groups: Optional[Dict[str, List[str]]] = None,
    top_n: Optional[int] = None,
) -> Dict[str, Any]:
    """Build deterministic, method-aware evidence before any AI call."""
    normalized_method = "ora" if method in {"hypergeometric", "ora"} else str(method).lower()
    if normalized_method not in _DEFAULT_EVIDENCE_LIMITS:
        raise ValueError(f"Unsupported AI evidence method: {method}")
    evidence_limit = _DEFAULT_EVIDENCE_LIMITS[normalized_method] if top_n is None else top_n
    if not isinstance(evidence_limit, int) or isinstance(evidence_limit, bool) or evidence_limit < 1:
        raise ValueError("AI evidence top_n must be a positive integer")
    evidence: Dict[str, Any] = {}
    databases: Dict[str, Any] = {}

    for database, frame in results.items():
        frame = frame if isinstance(frame, pd.DataFrame) else pd.DataFrame()
        database_name = str(database)
        database_info = {"status": "empty" if frame.empty else "available", "evidence_ids": []}
        selected: List[Dict[str, Any]] = []
        prefix = _database_prefix(database_name)
        tf_database = _is_tf_database(database_name, frame)

        if not frame.empty and normalized_method in {"gsea", "ora"} and not (normalized_method == "ora" and tf_database):
            if normalized_method == "gsea":
                ranked = _rank_rows(
                    frame,
                    ("padj", "FDR", "Adjusted_P_Value", "p.adjust"),
                    ("NES", "ES", "nes", "es"),
                )
                positive = [item for item in ranked if (_number(_first_value(item[1], ("NES", "nes", "ES", "es")), 0.0) or 0.0) >= 0]
                negative = [item for item in ranked if (_number(_first_value(item[1], ("NES", "nes", "ES", "es")), 0.0) or 0.0) < 0]
                selected_rows = positive[:evidence_limit] + negative[:evidence_limit]
                id_prefix = f"GSEA_{prefix}"
                for position, row in selected_rows:
                    nes = _number(_first_value(row, ("NES", "nes")))
                    gsea_values = {
                        "pval": _first_value(row, ("pval", "P_Value", "p_value")),
                        "padj": _first_value(row, ("padj", "FDR", "Adjusted_P_Value", "p.adjust")),
                        "ES": _first_value(row, ("ES", "es")),
                        "NES": nes,
                        "direction": "positive" if (nes or 0) >= 0 else "negative",
                        "size": _first_value(row, ("size", "setSize", "Gene_Count")),
                        "leading_edge_genes": _split_genes(_first_value(row, ("leadingEdge", "Lead_genes", "Leading_Edge", "Genes"))),
                    }
                    if tf_database:
                        gsea_values.update({
                            "target_count": _first_value(row, ("Target_Count", "Gene_Count", "size", "Count")),
                            "target_genes": _split_genes(_first_value(row, ("Target_Genes", "Genes", "geneID", "leadingEdge"))),
                            "source": _first_value(row, ("Library", "Source", "Database", "Context"), database_name),
                            "consistency_rank": _first_value(row, ("Consistency_Score", "Rank", "NES", "ES")),
                        })
                    selected.append(_row_record(
                        f"{id_prefix}:R{len(selected) + 1:03d}", database_name, normalized_method,
                        "tf" if tf_database else "gsea", position, row, gsea_values,
                    ))
            else:
                ranked = _rank_rows(
                    frame,
                    ("Adjusted_P_Value", "padj", "FDR", "p.adjust", "P_Value", "pvalue"),
                    ("EnrichFactor", "Rich_Factor", "RichFactor", "Gene_Ratio", "Gene_Count", "Count"),
                )[:evidence_limit]
                for position, row in ranked:
                    selected.append(_row_record(
                        f"{prefix}:R{len(selected) + 1:03d}", database_name, normalized_method, "ora",
                        position, row, {
                            "p_value": _first_value(row, ("P_Value", "pvalue", "p_value", "pval")),
                            "adjusted_p_value": _first_value(row, ("Adjusted_P_Value", "padj", "FDR", "p.adjust")),
                            "gene_count": _first_value(row, ("Gene_Count", "Count", "gene_count", "size")),
                            "genes": _split_genes(_first_value(row, ("Genes", "geneID", "gene_ids"))),
                            "enrich_factor": _first_value(row, ("EnrichFactor", "Rich_Factor", "RichFactor", "Gene_Ratio")),
                        },
                    ))
        elif not frame.empty and normalized_method in {"ssgsea", "gsva"}:
            sample_columns = _activity_sample_columns(frame)
            ranked = []
            for position, (_, row) in enumerate(frame.iterrows()):
                values = [_number(row.get(column)) for column in sample_columns]
                observed = [value for value in values if value is not None]
                activity_summary = _activity_summary(row, sample_columns, groups)
                differences = [
                    abs(_number(value, 0.0) or 0.0)
                    for value in activity_summary["group_differences"].values()
                ]
                if differences:
                    score = max(differences)
                else:
                    score = pd.Series(observed).var() if len(observed) > 1 else 0.0
                ranked.append((-score, position, row))
            ranked.sort(key=lambda item: (item[0], item[1]))
            id_prefix = f"{'ssGSEA' if normalized_method == 'ssgsea' else 'GSVA'}_{prefix}"
            for _, position, row in ranked[:evidence_limit]:
                selected.append(_row_record(
                    f"{id_prefix}:R{len(selected) + 1:03d}", database_name, normalized_method, "activity",
                    position, row, _activity_summary(row, sample_columns, groups),
                ))
        elif not frame.empty and tf_database:
            ranked = _rank_rows(
                frame,
                ("FDR", "padj", "Adjusted_P_Value", "p.adjust", "Pvalue", "P_Value", "pvalue"),
                ("NES", "Consistency_Score", "Rank", "Gene_Count", "Target_Count", "size"),
            )[:evidence_limit]
            id_prefix = f"TF_{prefix}"
            for position, row in ranked:
                selected.append(_row_record(
                    f"{id_prefix}:R{len(selected) + 1:03d}", database_name, normalized_method, "tf",
                    position, row, {
                        "p_value": _first_value(row, ("Pvalue", "P_Value", "pvalue", "pval")),
                        "fdr": _first_value(row, ("FDR", "padj", "Adjusted_P_Value", "p.adjust")),
                        "target_count": _first_value(row, ("Target_Count", "Gene_Count", "size", "Count")),
                        "target_genes": _split_genes(_first_value(row, ("Target_Genes", "Genes", "geneID", "leadingEdge"))),
                        "consistency_rank": _first_value(row, ("Consistency_Score", "Rank", "NES", "ES")),
                        "source": _first_value(row, ("Library", "Source", "Database", "Context"), database_name),
                    },
                ))

        for record in selected:
            evidence_id = record["evidence_id"]
            evidence[evidence_id] = record
            database_info["evidence_ids"].append(evidence_id)
        database_info["evidence_count"] = len(selected)
        databases[database_name] = database_info

    return {
        "schema_version": 1,
        "method": normalized_method,
        "databases": databases,
        "evidence": evidence,
    }


def build_interpretation_prompt(evidence: Dict[str, Any], mode: str = "summary", context: str = "") -> str:
    """Build the single shared prompt used by every AI backend."""
    if mode not in AI_MODES:
        raise ValueError(f"Unknown AI interpretation mode: {mode}. Available: {sorted(AI_MODES)}")
    mode_instruction = {
        "summary": "Write a concise evidence-based summary for a results report.",
        "reviewer": "Act as a critical reviewer; highlight over-interpretation and statistical limitations.",
        "caption": "Write concise text suitable for a paper figure caption.",
    }[mode]
    return f"""You are interpreting an AllEnricher enrichment analysis. {mode_instruction}
Use only the supplied evidence. Do not infer clinical meaning, causality, experimental procedures, or citations.
Return JSON only, with no Markdown fences, using exactly this top-level structure:
{{"schema_version":1,"method":"{evidence['method']}","profile":"{mode}","databases":{{}}}}
For every database include exactly these arrays: core_themes, key_evidence, limitations, computational_checks.
Each array item must be an object with a concise text field and an evidence_ids array.
Every data-based statement must cite one or more evidence IDs from the supplied evidence; never invent IDs.
Separate positive and negative GSEA directions when applicable. Mention leading-edge genes only when supplied.
Keep the response concise: use at most 3 items per section and keep each text under 35 words.
{context}

Structured evidence:
{json.dumps(evidence, ensure_ascii=False, indent=2)}"""


def _validate_section_items(database: str, section: str, items: Any, valid_ids: set[str]) -> List[Dict[str, Any]]:
    if not isinstance(items, list):
        raise ValueError(f"AI output field databases.{database}.{section} must be an array")
    validated = []
    for index, item in enumerate(items):
        if not isinstance(item, dict) or not isinstance(item.get("text"), str):
            raise ValueError(f"AI output field databases.{database}.{section}[{index}] must contain text")
        ids = item.get("evidence_ids")
        if not isinstance(ids, list) or not ids or not all(isinstance(value, str) for value in ids):
            raise ValueError(
                f"AI output field databases.{database}.{section}[{index}] must contain at least one evidence_id"
            )
        unknown = [value for value in ids if value not in valid_ids]
        if unknown:
            raise ValueError(f"AI output references unknown evidence_id(s): {', '.join(unknown)}")
        validated.append({"text": item["text"], "evidence_ids": ids})
    return validated


def validate_interpretation(payload: Any, evidence: Dict[str, Any], mode: str) -> Dict[str, Any]:
    """Validate model JSON and attach the code-generated evidence mapping."""
    if isinstance(payload, str):
        text = payload.strip()
        if text.startswith("```"):
            text = re.sub(r"^```(?:json)?\s*|\s*```$", "", text, flags=re.IGNORECASE | re.DOTALL).strip()
        try:
            payload = json.loads(text)
        except json.JSONDecodeError as exc:
            raise ValueError(f"AI output is not valid JSON: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError("AI output must be a JSON object")
    if payload.get("schema_version") != 1:
        raise ValueError("AI output schema_version must be 1")
    if payload.get("method") != evidence["method"]:
        raise ValueError("AI output method does not match the analysis method")
    if payload.get("profile") != mode:
        raise ValueError("AI output profile does not match --ai-mode")
    databases = payload.get("databases")
    expected = set(evidence["databases"])
    if not isinstance(databases, dict) or set(databases) != expected:
        missing = sorted(expected - set(databases or {}))
        extra = sorted(set(databases or {}) - expected)
        raise ValueError(f"AI output databases do not match evidence (missing={missing}, extra={extra})")
    valid_ids = set(evidence["evidence"])
    normalized_databases = {}
    for database in evidence["databases"]:
        source = databases[database]
        if not isinstance(source, dict):
            raise ValueError(f"AI output database '{database}' must be an object")
        normalized_databases[database] = {
            section: _validate_section_items(database, section, source.get(section), valid_ids)
            for section in _AI_SECTIONS
        }
    return {
        "schema_version": 1,
        "method": evidence["method"],
        "profile": mode,
        "databases": normalized_databases,
        "evidence": evidence["evidence"],
    }


class AIInterpreterBase(ABC):
    """Abstract interface for enrichment-result interpreters."""

    @abstractmethod
    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """Interpret result tables grouped by database."""
        pass

    @abstractmethod
    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """Summarize one biological term and its matched genes."""
        pass


class OpenAIInterpreter(AIInterpreterBase):
    """Interpret enrichment results with an OpenAI chat model."""

    def __init__(
        self,
        api_key: str = None,
        model: str = "gpt-4",
        max_tokens: int = 2000,
        temperature: float = 0.7
    ):
        """Configure the OpenAI client and generation settings."""
        self.api_key = api_key or os.getenv("OPENAI_API_KEY")
        self.model = model
        self.max_tokens = max_tokens
        self.temperature = temperature

        if not self.api_key:
            logger.warning("OpenAI API key not provided. AI interpretation will be disabled.")

    def _call_api(self, prompt: str) -> str:
        """Send one interpretation prompt to the OpenAI chat API."""
        try:
            import openai

            client = openai.OpenAI(api_key=self.api_key)

            response = client.chat.completions.create(
                model=self.model,
                messages=[
                    {
                        "role": "system",
                        "content": (
                            "You are a bioinformatics specialist interpreting gene set "
                            "enrichment tables. Use only the supplied results, distinguish "
                            "statistical evidence from biological hypotheses, and do not "
                            "invent pathways, genes, experimental details, or citations."
                        ),
                    },
                    {"role": "user", "content": prompt}
                ],
                max_tokens=self.max_tokens,
                temperature=self.temperature
            )

            return response.choices[0].message.content

        except ImportError:
            logger.error("openai package not installed. Run: pip install openai")
            return "Error: openai package not installed"
        except Exception as e:
            logger.error(f"OpenAI API error: {e}")
            return f"Error: {str(e)}"

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """Generate a concise interpretation for each non-empty database table."""
        interpretations = {}

        if not self.api_key:
            return interpretations

        for db_name, df in results.items():
            if len(df) == 0:
                interpretations[db_name] = "No enrichment terms were available for interpretation."
                continue

            top_results = df.head(20)
            summary_lines = []

            for _, row in top_results.iterrows():
                summary_lines.append(
                    f"- {row.get('Term_Name', 'N/A')}: "
                    f"P-value={row.get('P_Value', 1):.2e}, "
                    f"Genes={row.get('Gene_Count', 0)}"
                )

            prompt = f"""Interpret the following AllEnricher {db_name} result table in no more than 250 words.

Use only the terms, P values, and gene counts shown below. Do not infer the
experimental design or claim causality. Describe broad patterns as hypotheses
that require domain review and literature validation.

Top {len(top_results)} enriched terms:
{chr(10).join(summary_lines)}

Use exactly these headings:
**Main themes**: [one or two evidence-based sentences]
**Representative terms**: [comma-separated term names]
**Statistical context**: [one sentence limited to the supplied statistics]
**Gene-level observations**: [one sentence, or state that gene details are insufficient]"""

            interpretation = self._call_api(prompt)
            interpretations[db_name] = interpretation

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """Generate a short description of one enriched term."""
        if not self.api_key:
            return ""

        prompt = f"""
Briefly describe the established biological meaning of this database term.
Do not infer the experimental context or claim that the matched genes prove a
mechanism.

Term: {term_name}
Associated genes: {', '.join(gene_list[:10])}{'...' if len(gene_list) > 10 else ''}

Write two or three concise sentences and use only the supplied term and genes.
"""

        return self._call_api(prompt)


class ClaudeInterpreter(AIInterpreterBase):
    """Interpret enrichment results with an Anthropic Claude model."""

    def __init__(
        self,
        api_key: str = None,
        model: str = "claude-3-opus-20240229",
        max_tokens: int = 2000
    ):
        """Configure the Anthropic client and model."""
        self.api_key = api_key or os.getenv("ANTHROPIC_API_KEY")
        self.model = model
        self.max_tokens = max_tokens

        if not self.api_key:
            logger.warning("Anthropic API key not provided. AI interpretation will be disabled.")

    def _call_api(self, prompt: str) -> str:
        """Send one interpretation prompt to the Anthropic Messages API."""
        try:
            import anthropic

            client = anthropic.Anthropic(api_key=self.api_key)

            message = client.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                messages=[
                    {"role": "user", "content": prompt}
                ]
            )

            return message.content[0].text

        except ImportError:
            logger.error("anthropic package not installed. Run: pip install anthropic")
            return "Error: anthropic package not installed"
        except Exception as e:
            logger.error(f"Anthropic API error: {e}")
            return f"Error: {str(e)}"

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """Generate a concise interpretation for each non-empty database table."""
        interpretations = {}

        if not self.api_key:
            return interpretations

        for db_name, df in results.items():
            if len(df) == 0:
                interpretations[db_name] = "No enrichment terms were available for interpretation."
                continue

            top_results = df.head(20)
            summary_lines = []

            for _, row in top_results.iterrows():
                summary_lines.append(
                    f"- {row.get('Term_Name', 'N/A')}: "
                    f"P-value={row.get('P_Value', 1):.2e}"
                )

            prompt = f"""Interpret the following AllEnricher {db_name} result table in no more than 300 words.

Use only the supplied terms and P values. Do not infer an experimental design,
claim causality, or invent genes, pathways, or citations.

Top {len(top_results)} enriched terms:
{chr(10).join(summary_lines)}

Use exactly these headings:
**Main themes**
**Representative terms**
**Statistical context**
**Limitations**

End by stating that the interpretation requires domain review and literature validation."""

            interpretation = self._call_api(prompt)
            interpretations[db_name] = interpretation

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """Generate a short description of one enriched term."""
        if not self.api_key:
            return ""

        prompt = (
            f"Describe the established biological meaning of the database term '{term_name}' "
            "in two or three concise sentences. Do not infer experimental context or causality."
        )
        return self._call_api(prompt)


class OllamaInterpreter(AIInterpreterBase):
    """Interpret enrichment results with a locally hosted Ollama model."""

    def __init__(self, model: str = "llama2", base_url: str = "http://localhost:11434"):
        """Configure the Ollama model and service URL."""
        self.model = model
        self.base_url = base_url

    def _call_api(self, prompt: str) -> str:
        """Send one interpretation prompt to the local Ollama service."""
        try:
            import requests

            response = requests.post(
                f"{self.base_url}/api/generate",
                json={
                    "model": self.model,
                    "prompt": prompt,
                    "stream": False,
                }
            )

            if response.status_code == 200:
                return response.json().get("response", "")
            else:
                return f"Error: {response.status_code}"

        except ImportError:
            logger.error("requests package not installed. Run: pip install requests")
            return "Error: requests package not installed"
        except Exception as e:
            logger.error(f"Ollama API error: {e}")
            return f"Error: {str(e)}"

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """Generate a concise interpretation for each non-empty database table."""
        interpretations = {}

        for db_name, df in results.items():
            if len(df) == 0:
                interpretations[db_name] = "No enrichment terms were available for interpretation."
                continue

            top_results = df.head(20)
            summary_lines = []

            for _, row in top_results.iterrows():
                summary_lines.append(f"- {row.get('Term_Name', 'N/A')}")

            prompt = f"""Interpret the following AllEnricher {db_name} result table in no more than 250 words.

Only term names are provided. Do not infer statistical strength, experimental
design, causality, genes, or citations that are not shown.

Top {len(top_results)} enriched terms:
{chr(10).join(summary_lines)}

Use exactly these headings:
**Main themes**
**Representative terms**
**Limitations**"""

            interpretation = self._call_api(prompt)
            interpretations[db_name] = interpretation

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """Generate a short local-model description of one term."""
        prompt = (
            f"Describe the established biological meaning of the database term '{term_name}' "
            "in two concise sentences without inferring experimental context."
        )
        return self._call_api(prompt)


class DeepSeekInterpreter(AIInterpreterBase):
    """Interpret enrichment results with the DeepSeek API."""

    # DeepSeek API Base URL
    BASE_URL = "https://api.deepseek.com"

    def __init__(
        self,
        api_key: str = None,
        model: str = "deepseek-chat",
        max_tokens: int = 2000,
        temperature: float = 0.7
    ):
        """Configure the DeepSeek client and generation settings."""
        self.api_key = api_key or os.getenv("DEEPSEEK_API_KEY")
        self.model = model
        self.max_tokens = max_tokens
        self.temperature = temperature

        if not self.api_key:
            logger.warning("DeepSeek API key not provided. Set DEEPSEEK_API_KEY environment variable.")

    def _call_api(self, prompt: str) -> str:
        """Send one interpretation prompt to the DeepSeek chat API."""
        try:
            import openai

            # DeepSeek exposes an OpenAI-compatible chat endpoint.
            client = openai.OpenAI(
                api_key=self.api_key,
                base_url=self.BASE_URL
            )

            response = client.chat.completions.create(
                model=self.model,
                messages=[
                    {
                        "role": "system",
                        "content": (
                            "You are a bioinformatics analyst interpreting gene-set enrichment "
                            "results. Base every statement on the supplied terms and statistics. "
                            "Do not invent genes, pathways, experimental conditions, causal "
                            "claims, or references. Clearly identify limitations in the evidence."
                        ),
                    },
                    {"role": "user", "content": prompt}
                ],
                max_tokens=self.max_tokens,
                temperature=self.temperature
            )

            return response.choices[0].message.content

        except ImportError:
            logger.error("openai package not installed. Run: pip install openai")
            return "Error: openai package not installed"
        except Exception as e:
            logger.error(f"DeepSeek API error: {e}")
            return f"Error: {str(e)}"

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """Generate interpretations using the shared result-table workflow."""
        interpretations = {}

        if not self.api_key:
            return interpretations

        for db_name, df in results.items():
            if len(df) == 0:
                interpretations[db_name] = "No enrichment terms were available for interpretation."
                continue

            top_results = df.head(20)
            summary_lines = []

            for _, row in top_results.iterrows():
                summary_lines.append(
                    f"- {row.get('Term_Name', 'N/A')}: "
                    f"P value={row.get('P_Value', 1):.2e}, "
                    f"gene count={row.get('Gene_Count', 0)}"
                )

            prompt = f"""Interpret the following {db_name} enrichment results in no more than 250 words.

Top {len(top_results)} enriched terms:
{chr(10).join(summary_lines)}

Use exactly these headings:
**Main themes**: Summarize the dominant biological themes in one or two sentences.
**Key terms**: List only terms present above.
**Interpretation**: Explain what the enrichment supports without inferring causality.
**Limitations**: State what cannot be concluded from these data alone."""

            interpretations[db_name] = self._call_api(prompt)

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """Generate a short description of one enriched term."""
        if not self.api_key:
            return ""

        genes = ", ".join(gene_list[:10]) or "not provided"
        prompt = f"""Describe the biological meaning of the enriched term "{term_name}" in two or three sentences.
Input genes associated with this term: {genes}.
Use only the supplied term and genes. Do not infer the experimental context or cite unverified sources."""

        return self._call_api(prompt)


class GLMInterpreter(AIInterpreterBase):
    """Interpret enrichment results with a GLM-compatible API."""

    # OpenAI-compatible GLM endpoint.
    BASE_URL = "https://open.bigmodel.cn/api/paas/v4"

    def __init__(
        self,
        api_key: str = None,
        model: str = "glm-4",
        max_tokens: int = 2000,
        temperature: float = 0.7
    ):
        """Configure the GLM client and generation settings."""
        self.api_key = api_key or os.getenv("GLM_API_KEY")
        self.model = model
        self.max_tokens = max_tokens
        self.temperature = temperature

        if not self.api_key:
            logger.warning("GLM API key not provided. Set GLM_API_KEY environment variable.")

    def _call_api(self, prompt: str) -> str:
        """Send one interpretation prompt to the GLM chat API."""
        try:
            import openai

            client = openai.OpenAI(
                api_key=self.api_key,
                base_url=self.BASE_URL
            )

            response = client.chat.completions.create(
                model=self.model,
                messages=[
                    {
                        "role": "system",
                        "content": (
                            "You are a bioinformatics analyst interpreting gene-set enrichment "
                            "results. Use only the supplied terms and statistics. Do not invent "
                            "genes, experimental conditions, causal claims, or references."
                        ),
                    },
                    {"role": "user", "content": prompt}
                ],
                max_tokens=self.max_tokens,
                temperature=self.temperature
            )

            return response.choices[0].message.content

        except ImportError:
            logger.error("openai package not installed. Run: pip install openai")
            return "Error: openai package not installed"
        except Exception as e:
            logger.error(f"GLM API error: {e}")
            return f"Error: {str(e)}"

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """Generate interpretations using the shared result-table workflow."""
        interpretations = {}

        if not self.api_key:
            return interpretations

        for db_name, df in results.items():
            if len(df) == 0:
                interpretations[db_name] = "No enrichment terms were available for interpretation."
                continue

            top_results = df.head(20)
            summary_lines = []

            for _, row in top_results.iterrows():
                summary_lines.append(
                    f"- {row.get('Term_Name', 'N/A')}: "
                    f"P value={row.get('P_Value', 1):.2e}, "
                    f"gene count={row.get('Gene_Count', 0)}"
                )

            prompt = f"""Interpret the following {db_name} enrichment results in no more than 250 words.

Top {len(top_results)} enriched terms:
{chr(10).join(summary_lines)}

Use exactly these headings:
**Main themes**: Summarize the dominant biological themes in one or two sentences.
**Key terms**: List only terms present above.
**Interpretation**: Explain what the enrichment supports without inferring causality.
**Limitations**: State what cannot be concluded from these data alone."""

            interpretations[db_name] = self._call_api(prompt)

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """Generate a short description of one enriched term."""
        if not self.api_key:
            return ""

        genes = ", ".join(gene_list[:10]) or "not provided"
        prompt = f"""Describe the biological meaning of the enriched term "{term_name}" in two or three sentences.
Input genes associated with this term: {genes}.
Use only the supplied term and genes. Do not infer the experimental context or cite unverified sources."""
        return self._call_api(prompt)


class MiniMaxInterpreter(AIInterpreterBase):
    """Interpret enrichment results with the MiniMax API."""

    # OpenAI-compatible MiniMax endpoint.
    BASE_URL = "https://api.minimax.chat/v1"

    def __init__(
        self,
        api_key: str = None,
        group_id: str = None,
        model: str = "abab6.5s-chat",
        max_tokens: int = 2000,
        temperature: float = 0.7
    ):
        """Configure the MiniMax client, model, and group identifier."""
        self.api_key = api_key or os.getenv("MINIMAX_API_KEY")
        self.group_id = group_id or os.getenv("MINIMAX_GROUP_ID")
        self.model = model
        self.max_tokens = max_tokens
        self.temperature = temperature

        if not self.api_key:
            logger.warning("MiniMax API key not provided. Set MINIMAX_API_KEY environment variable.")
        if not self.group_id:
            logger.warning("MiniMax Group ID not provided. Set MINIMAX_GROUP_ID environment variable.")

    def _call_api(self, prompt: str) -> str:
        """Send one interpretation prompt to the MiniMax chat API."""
        try:
            import openai

            if not self.api_key or not self.group_id:
                return "Error: MiniMax API key or Group ID not configured"

            client = openai.OpenAI(
                api_key=self.api_key,
                base_url=f"{self.BASE_URL}"
            )

            response = client.chat.completions.create(
                model=self.model,
                messages=[
                    {
                        "role": "system",
                        "content": (
                            "You are a bioinformatics analyst interpreting gene-set enrichment "
                            "results. Use only the supplied terms and statistics. Do not invent "
                            "genes, experimental conditions, causal claims, or references."
                        ),
                    },
                    {"role": "user", "content": prompt}
                ],
                max_tokens=self.max_tokens,
                temperature=self.temperature
            )

            return response.choices[0].message.content

        except ImportError:
            logger.error("openai package not installed. Run: pip install openai")
            return "Error: openai package not installed"
        except Exception as e:
            logger.error(f"MiniMax API error: {e}")
            return f"Error: {str(e)}"

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """Generate interpretations using the shared result-table workflow."""
        interpretations = {}

        if not self.api_key or not self.group_id:
            return interpretations

        for db_name, df in results.items():
            if len(df) == 0:
                interpretations[db_name] = "No enrichment terms were available for interpretation."
                continue

            top_results = df.head(20)
            summary_lines = []

            for _, row in top_results.iterrows():
                summary_lines.append(
                    f"- {row.get('Term_Name', 'N/A')}: "
                    f"P-value={row.get('P_Value', 1):.2e}, "
                    f"Genes={row.get('Gene_Count', 0)}"
                )

            prompt = f"""Interpret the following {db_name} enrichment results in no more than 250 words.

Top {len(top_results)} enriched terms:
{chr(10).join(summary_lines)}

Use exactly these headings:
**Main themes**: Summarize the dominant biological themes in one or two sentences.
**Key terms**: List only terms present above.
**Interpretation**: Explain what the enrichment supports without inferring causality.
**Limitations**: State what cannot be concluded from these data alone."""

            interpretations[db_name] = self._call_api(prompt)

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """Generate a short description of one enriched term."""
        if not self.api_key or not self.group_id:
            return ""

        genes = ", ".join(gene_list[:10]) or "not provided"
        prompt = f"""Describe the biological meaning of the enriched term "{term_name}" in two or three sentences.
Input genes associated with this term: {genes}.
Use only the supplied term and genes. Do not infer the experimental context or cite unverified sources."""
        return self._call_api(prompt)


class MockInterpreter(AIInterpreterBase):
    """Deterministic interpreter used by tests and offline examples."""

    def interpret(self, results: Dict[str, pd.DataFrame], context: str = "") -> Dict[str, str]:
        """Return deterministic summaries without calling an external service."""
        interpretations = {}

        for db_name, df in results.items():
            if len(df) == 0:
                continue

            # Keep the deterministic fixture compact while preserving term names.
            top_terms = df.head(20)['Term_Name'].tolist()

            interpretation = f"""Offline test interpretation for {db_name} ({len(df)} terms analyzed; up to 20 listed)

**Main themes**: This deterministic test backend does not infer enrichment themes.

**Key terms**: See the enriched terms listed below.

**Interpretation**: No biological interpretation was generated because no external language model was used.

**Limitations**: This text validates report rendering only and must not be used as a scientific interpretation.

Top {len(top_terms)} enriched terms:
{chr(10).join([f"- {term}" for term in top_terms])}"""

            interpretations[db_name] = interpretation

        return interpretations

    def summarize_term(self, term_name: str, gene_list: List[str]) -> str:
        """Return a deterministic term summary for tests."""
        return f"The term '{term_name}' is associated with {len(gene_list)} genes from your input set."

    def structured_response(self, evidence: Dict[str, Any], mode: str) -> Dict[str, Any]:
        """Return a deterministic response with the same contract as real backends."""
        databases = {}
        for database, info in evidence["databases"].items():
            ids = list(info["evidence_ids"])
            first_ids = ids[:2]
            databases[database] = {
                "core_themes": ([{
                    "text": f"The {evidence['method']} evidence for {database} is represented by the selected records.",
                    "evidence_ids": first_ids,
                }] if first_ids else []),
                "key_evidence": [
                    {"text": f"Representative record {evidence_id} is available for review.", "evidence_ids": [evidence_id]}
                    for evidence_id in ids[:3]
                ],
                "limitations": ([{
                    "text": "This deterministic response does not establish biological causality.",
                    "evidence_ids": [ids[0]],
                }] if ids else []),
                "computational_checks": ([{
                    "text": "Review the linked source row and test threshold stability before interpretation.",
                    "evidence_ids": [ids[0]],
                }] if ids else []),
            }
        return {
            "schema_version": 1,
            "method": evidence["method"],
            "profile": mode,
            "databases": databases,
        }


class AIInterpreter:
    """Facade that constructs and delegates to the selected backend."""

    # Public backend name to implementation class.
    BACKENDS = {
        "openai": OpenAIInterpreter,
        "claude": ClaudeInterpreter,
        "deepseek": DeepSeekInterpreter,
        "glm": GLMInterpreter,
        "minimax": MiniMaxInterpreter,
        "ollama": OllamaInterpreter,
        "mock": MockInterpreter
    }

    def __init__(
        self,
        backend: str = "openai",
        api_key: str = None,
        model: str = None,
        **kwargs
    ):
        """Initialize one supported interpretation backend."""
        if backend not in self.BACKENDS:
            raise ValueError(f"Unknown backend: {backend}. Available: {list(self.BACKENDS.keys())}")

        self.backend_name = backend
        interpreter_class = self.BACKENDS[backend]

        # Each provider has a small set of backend-specific constructor options.
        if backend == "openai":
            self.interpreter = interpreter_class(
                api_key=api_key,
                model=model or "gpt-4",
                **kwargs
            )
        elif backend == "claude":
            self.interpreter = interpreter_class(
                api_key=api_key,
                model=model or "claude-3-opus-20240229",
                **kwargs
            )
        elif backend == "deepseek":
            self.interpreter = interpreter_class(
                api_key=api_key,
                model=model or "deepseek-chat",
                **kwargs
            )
        elif backend == "glm":
            self.interpreter = interpreter_class(
                api_key=api_key,
                model=model or "glm-4",
                **kwargs
            )
        elif backend == "minimax":
            # MiniMax additionally requires a group identifier.
            group_id = kwargs.pop("group_id", None)
            self.interpreter = interpreter_class(
                api_key=api_key,
                group_id=group_id,
                model=model or "abab6.5s-chat",
                **kwargs
            )
        elif backend == "ollama":
            self.interpreter = interpreter_class(
                model=model or "llama2",
                **kwargs
            )
        else:
            # The deterministic mock backend requires no credentials.
            self.interpreter = interpreter_class(**kwargs)
    
    def interpret_results(
        self,
        results: Dict[str, pd.DataFrame],
        context: str = "",
        include_term_summaries: bool = False
    ) -> Dict[str, Any]:
        """Interpret result tables and return an empty mapping on backend failure."""
        interpretations = self.interpreter.interpret(results, context)

        # Optional per-term summaries are limited to the first 20 result rows.
        if include_term_summaries:
            for db_name, df in results.items():
                if len(df) > 0:
                    term_summaries = {}
                    for _, row in df.head(20).iterrows():
                        term_name = row.get('Term_Name', '')
                        # The canonical result schema stores genes as semicolon-separated IDs.
                        genes = row.get('Genes', '').split(';')
                        if term_name:
                            term_summaries[term_name] = self.interpreter.summarize_term(term_name, genes)

                    interpretations[f"{db_name}_term_summaries"] = term_summaries

        return interpretations

    def interpret_structured_results(
        self,
        results: Dict[str, pd.DataFrame],
        method: str = "hypergeometric",
        mode: str = "summary",
        groups: Optional[Dict[str, List[str]]] = None,
        context: str = "",
        top_n: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Build evidence, call one backend, and validate its structured JSON."""
        if mode not in AI_MODES:
            raise ValueError(f"Unknown AI interpretation mode: {mode}. Available: {sorted(AI_MODES)}")
        evidence = build_structured_evidence(results, method=method, groups=groups, top_n=top_n)
        if self.backend_name == "mock":
            payload = self.interpreter.structured_response(evidence, mode)
        else:
            prompt = build_interpretation_prompt(evidence, mode=mode, context=context)
            raw = self.interpreter._call_api(prompt)
            if not isinstance(raw, str) or raw.startswith("Error:"):
                raise ValueError(f"AI backend failed: {raw}")
            payload = raw
        return validate_interpretation(payload, evidence, mode)
    
    def generate_report_section(
        self,
        results: Dict[str, pd.DataFrame],
        context: str = ""
    ) -> str:
        """Render interpretation text as an HTML report section."""
        interpretations = self.interpret_results(results, context)

        # Begin with an explicit review disclaimer.
        html_parts = ['''
        <div class="section" id="ai-interpretation">
            <h2><i class="fas fa-brain"></i> AI-Powered Interpretation</h2>
            <p class="ai-disclaimer">
                <i class="fas fa-info-circle"></i>
                The following interpretations are generated by AI ({}) and should be reviewed by domain experts.
            </p>
        '''.format(self.backend_name)]

        for db_name, interpretation in interpretations.items():
            # Per-term summaries are data for other clients, not report sections.
            if db_name.endswith('_term_summaries'):
                continue

            html_parts.append(f'''
            <div class="ai-interpretation">
                <h3><i class="fas fa-robot"></i> {db_name}</h3>
                <div class="interpretation-content">
                    {interpretation.replace(chr(10), '<br>')}
                </div>
            </div>
            ''')

        html_parts.append('</div>')
        return ''.join(html_parts)


def get_available_backends() -> List[str]:
    """Return the supported interpreter backend names."""
    return list(AIInterpreter.BACKENDS.keys())


def create_interpreter(
    backend: str = "mock",
    api_key: str = None,
    model: str = None,
    **kwargs
) -> AIInterpreter:
    """Create the facade for one backend.

    Backend-specific keyword arguments are forwarded without changing analysis
    results or configuration outside the interpreter."""
    return AIInterpreter(backend=backend, api_key=api_key, model=model, **kwargs)


def create_interpreter_from_config(config, backend: str = None) -> AIInterpreter:
    """Create an interpreter from an AllEnricher configuration object."""
    backend_name = backend or getattr(config, 'ai_backend', 'mock')

    # Backend-specific configuration takes precedence over global settings.
    api_key = config.get_ai_api_key(backend_name)

    # Resolve the model after the backend so provider defaults remain explicit.
    default_models = {
        "openai": "gpt-4",
        "claude": "claude-3-opus-20240229",
        "deepseek": "deepseek-chat",
        "glm": "glm-4",
        "minimax": "abab6.5s-chat",
        "ollama": "llama2",
        "mock": None,
    }
    model = config.get_ai_model(backend_name, default_models.get(backend_name))

    kwargs = {}

    # Custom endpoint, primarily used by Ollama-compatible deployments.
    base_url = config.get_ai_base_url(backend_name)
    if base_url:
        kwargs['base_url'] = base_url

    # MiniMax provider account identifier.
    group_id = config.get_ai_group_id(backend_name)
    if group_id:
        kwargs['group_id'] = group_id

    return AIInterpreter(backend=backend_name, api_key=api_key, model=model, **kwargs)
