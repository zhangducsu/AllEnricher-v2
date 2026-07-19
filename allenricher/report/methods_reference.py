"""Build publication-ready Methods text from recorded analysis metadata."""

from __future__ import annotations

from html import escape
from typing import Any, Dict, Iterable, List, Mapping, Optional

from allenricher.core.config import database_catalog_entry, database_display_name


ALLENRICHER_REFERENCE = {
    "source": "AllEnricher",
    "citation": (
        "Zhang D, Hu Q, Liu X, Zou K, Sarkodie EK, Liu X, et al. "
        "AllEnricher: a comprehensive gene set function enrichment tool for both "
        "model and non-model species. BMC Bioinformatics. 2020;21:106. "
        "doi:10.1186/s12859-020-3408-y."
    ),
    "url": "https://doi.org/10.1186/s12859-020-3408-y",
}


DATABASE_REFERENCES: Dict[str, List[Dict[str, str]]] = {
    "GO": [
        {
            "source": "GO",
            "citation": (
                "Ashburner M, Ball CA, Blake JA, et al. Gene ontology: tool for the "
                "unification of biology. Nat Genet. 2000;25(1):25-29. "
                "doi:10.1038/75556."
            ),
            "url": "https://doi.org/10.1038/75556",
        },
        {
            "source": "GO",
            "citation": (
                "The Gene Ontology Consortium. The Gene Ontology knowledgebase in "
                "2026. Nucleic Acids Res. 2026;54(D1):D1779-D1792. "
                "doi:10.1093/nar/gkaf1292."
            ),
            "url": "https://doi.org/10.1093/nar/gkaf1292",
        },
    ],
    "KEGG": [
        {
            "source": "KEGG",
            "citation": (
                "Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. "
                "Nucleic Acids Res. 2000;28(1):27-30. doi:10.1093/nar/28.1.27."
            ),
            "url": "https://doi.org/10.1093/nar/28.1.27",
        }
    ],
    "REACTOME": [
        {
            "source": "Reactome",
            "citation": (
                "Ragueneau E, Gong C, Sinquin P, et al. The Reactome Knowledgebase "
                "2026. Nucleic Acids Res. 2025 Nov 18. "
                "doi:10.1093/nar/gkaf1223."
            ),
            "url": "https://doi.org/10.1093/nar/gkaf1223",
        }
    ],
    "WIKIPATHWAYS": [
        {
            "source": "WikiPathways",
            "citation": (
                "Agrawal A, Balci H, Hanspers K, et al. WikiPathways 2024: next "
                "generation pathway database. Nucleic Acids Res. "
                "2024;52(D1):D679-D689. doi:10.1093/nar/gkad960."
            ),
            "url": "https://doi.org/10.1093/nar/gkad960",
        }
    ],
    "DO": [
        {
            "source": "Disease Ontology",
            "citation": (
                "Baron JA, Sanchez-Beato Johnson CM, Schor MA, et al. The DO-KB "
                "knowledgebase 2026 update: expanding programmatic and language "
                "access. Nucleic Acids Res. 2026;54(D1):D1425-D1436. "
                "doi:10.1093/nar/gkaf1213."
            ),
            "url": "https://doi.org/10.1093/nar/gkaf1213",
        }
    ],
    "DISGENET": [
        {
            "source": "DisGeNET",
            "citation": (
                "Pinero J, Ramirez-Anguita JM, Sauch-Pitarch J, et al. The DisGeNET "
                "knowledge platform for disease genomics: 2019 update. Nucleic "
                "Acids Res. 2020;48(D1):D845-D855. doi:10.1093/nar/gkz1021."
            ),
            "url": "https://doi.org/10.1093/nar/gkz1021",
        }
    ],
    "TRRUST": [
        {
            "source": "TRRUST",
            "citation": (
                "Han H, Cho JW, Lee S, et al. TRRUST v2: an expanded reference "
                "database of human and mouse transcriptional regulatory interactions. "
                "Nucleic Acids Res. 2018;46(D1):D380-D386. "
                "doi:10.1093/nar/gkx1013."
            ),
            "url": "https://doi.org/10.1093/nar/gkx1013",
        }
    ],
    "CHEA3": [
        {
            "source": "ChEA3",
            "citation": (
                "Keenan AB, Torre D, Lachmann A, et al. ChEA3: transcription factor "
                "enrichment analysis by orthogonal omics integration. Nucleic Acids "
                "Res. 2019;47(W1):W212-W224. doi:10.1093/nar/gkz446."
            ),
            "url": "https://doi.org/10.1093/nar/gkz446",
        }
    ],
    "ANIMALTFDB": [
        {
            "source": "AnimalTFDB",
            "citation": (
                "Shen WK, Chen SY, Gan ZQ, et al. AnimalTFDB 4.0: a comprehensive "
                "animal transcription factor database updated with variation and "
                "expression annotations. Nucleic Acids Res. 2023;51(D1):D39-D45. "
                "doi:10.1093/nar/gkac907."
            ),
            "url": "https://doi.org/10.1093/nar/gkac907",
        }
    ],
    "HTFTARGET": [
        {
            "source": "hTFtarget",
            "citation": (
                "Zhang Q, Liu W, Zhang HM, et al. hTFtarget: A Comprehensive Database "
                "for Regulations of Human Transcription Factors and Their Targets. "
                "Genomics Proteomics Bioinformatics. 2020;18(2):120-128. "
                "doi:10.1016/j.gpb.2019.09.006."
            ),
            "url": "https://doi.org/10.1016/j.gpb.2019.09.006",
        }
    ],
}


METHOD_NAMES = {
    "hypergeometric": "over-representation analysis (ORA)",
    "gsea": "gene set enrichment analysis (GSEA)",
    "ssgsea": "single-sample gene set enrichment analysis (ssGSEA)",
    "gsva": "gene set variation analysis (GSVA)",
}


def _language(value: Optional[str]) -> str:
    """Return the only supported report language.

    The argument remains accepted for API and CLI backward compatibility.
    """

    return "en"


def _display_value(value: Any, language: str) -> str:
    if value is None or value == "":
        return "To be added"
    if (isinstance(value, float) and value == float("inf")) or str(value).lower() in {
        "inf",
        "infinity",
        "unbounded",
    }:
        return "unbounded"
    return str(value)


def _species_text(metadata: Mapping[str, Any], language: str) -> str:
    code = _display_value(metadata.get("species"), language)
    name = metadata.get("species_name")
    taxid = metadata.get("species_taxonomy_id")
    if name and taxid:
        return f"{name} (TaxID: {taxid}; {code})"
    if name:
        return f"{name} ({code})"
    if taxid:
        return f"TaxID: {taxid} ({code})"
    return code


def _case_insensitive_get(values: Mapping[str, Any], key: str) -> Any:
    target = key.upper()
    for candidate, value in values.items():
        if str(candidate).upper() == target:
            return value
    return None


def _database_texts(metadata: Mapping[str, Any], language: str) -> List[str]:
    databases = list(metadata.get("databases") or [])
    versions = metadata.get("database_versions") or {}
    source_versions = metadata.get("source_versions") or {}
    items = []
    for database in databases:
        catalog = database_catalog_entry(str(database))
        version = (
            _case_insensitive_get(source_versions, str(database))
            or catalog.get("source_version")
            or _case_insensitive_get(versions, str(database))
        )
        version_text = _display_value(version, language)
        label = database_display_name(str(database))
        if catalog.get("source_version") == version and label != str(database):
            items.append(label)
            continue
        items.append(f"{label} (version: {version_text})")
    return items


def _join_items(items: Iterable[str], language: str) -> str:
    values = list(items)
    if not values:
        return "To be added"
    if len(values) > 1:
        return ", ".join(values[:-1]) + f", and {values[-1]}"
    return ", ".join(values)


def _size_text(parameters: Mapping[str, Any], language: str) -> str:
    limits = parameters.get("gene_set_size_by_database") or {}
    if not limits:
        return ""
    parts = []
    for database, values in limits.items():
        minimum = _display_value((values or {}).get("min"), language)
        maximum = _display_value((values or {}).get("max"), language)
        parts.append(f"{database_display_name(str(database))}: {minimum}-{maximum}")
    return " Gene-set size limits were " + "; ".join(parts) + "."


def _parameter_text(method: str, parameters: Mapping[str, Any], language: str) -> str:
    size_text = _size_text(parameters, language)
    if method == "hypergeometric":
        background = _display_value(parameters.get("background_mode"), language)
        correction = _display_value(parameters.get("correction"), language)
        pvalue = _display_value(parameters.get("pvalue_cutoff"), language)
        qvalue = _display_value(parameters.get("qvalue_cutoff"), language)
        return (
            f" The recorded background mode was {background}; multiple-testing "
            f"correction used {correction}, with P-value and adjusted P-value "
            f"cutoffs of {pvalue} and {qvalue}, respectively.{size_text}"
        )
    if method == "gsea":
        pvalue = _display_value(parameters.get("pvalue_cutoff"), language)
        qvalue = _display_value(parameters.get("qvalue_cutoff"), language)
        return f" P-value and adjusted P-value cutoffs were {pvalue} and {qvalue}, respectively.{size_text}"
    if method == "gsva":
        gsva_method = _display_value(parameters.get("gsva_method"), language)
        kcdf = _display_value(parameters.get("gsva_kcdf"), language)
        tau = _display_value(parameters.get("gsva_tau"), language)
        return f" The recorded GSVA settings were method={gsva_method}, kcdf={kcdf}, and tau={tau}.{size_text}"
    if method == "ssgsea":
        tau = _display_value(parameters.get("ssgsea_tau"), language)
        return f" The recorded ssGSEA tau parameter was {tau}.{size_text}"
    return size_text


def _reference_entries(databases: Iterable[str], language: str) -> List[Dict[str, str]]:
    references = [dict(ALLENRICHER_REFERENCE)]
    seen = {ALLENRICHER_REFERENCE["citation"]}
    for database in databases:
        entries = DATABASE_REFERENCES.get(str(database).upper())
        if not entries:
            missing = {
                "source": str(database),
                "citation": f"{database}: To be added",
                "url": "",
            }
            entries = [missing]
        for entry in entries:
            if entry["citation"] not in seen:
                references.append(dict(entry))
                seen.add(entry["citation"])
    return references


def build_methods_reference(
    metadata: Optional[Mapping[str, Any]],
    output_language: Optional[str] = None,
) -> Dict[str, Any]:
    """Return Methods prose and verified references from recorded metadata only."""

    metadata = metadata or {}
    language = _language(output_language or metadata.get("methods_language"))
    method = str(metadata.get("analysis_method") or "").lower()
    method_name = METHOD_NAMES.get(method)
    if not method_name:
        method_name = "To be added"
    version = _display_value(metadata.get("allenricher_version"), language)
    species = _species_text(metadata, language)
    databases = list(metadata.get("databases") or [])
    database_text = _join_items(_database_texts(metadata, language), language)
    parameters = metadata.get("parameters") or {}

    paragraph = (
        f"AllEnricher version {version} was used to perform {method_name} for "
        f"{species}. The analysis used {database_text}."
        f"{_parameter_text(method, parameters, language)}"
    )
    labels = {
        "title": "Materials and Methods Writing Reference",
        "analysis_heading": "Analysis Methods",
        "references_heading": "References",
    }

    return {
        **labels,
        "language": language,
        "paragraphs": [paragraph.strip()],
        "references": _reference_entries(databases, language),
    }


def render_methods_reference_html(
    metadata: Optional[Mapping[str, Any]],
    output_language: Optional[str] = None,
) -> str:
    """Render the shared Methods reference as a report section."""

    content = build_methods_reference(metadata, output_language)
    paragraphs = "".join(f"<p>{escape(text)}</p>" for text in content["paragraphs"])
    references = []
    for item in content["references"]:
        citation = escape(item["citation"])
        if item.get("url"):
            citation += (
                f' <a href="{escape(item["url"], quote=True)}" target="_blank" '
                'rel="noopener noreferrer">DOI</a>'
            )
        references.append(f"<li>{citation}</li>")
    return (
        '<section id="methods-reference" class="section methods-reference">'
        f'<h2>{escape(content["title"])}</h2>'
        f'<h3>{escape(content["analysis_heading"])}</h3>'
        f'<div class="methods-prose">{paragraphs}</div>'
        f'<h3>{escape(content["references_heading"])}</h3>'
        f'<ol class="methods-references">{"".join(references)}</ol>'
        '</section>'
    )
