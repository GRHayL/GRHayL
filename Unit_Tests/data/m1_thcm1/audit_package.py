#!/usr/bin/env python3
"""Audit the repository-local M1 fixture package without external inputs."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re


MANIFEST_NAME = "package_manifest.json"
SCHEMA = "grhayl_m1_thcm1_fixture_package_manifest"
VERSION = 1
PAYLOAD_SUFFIXES = frozenset((".m1", ".dat"))
STORED_REFERENCE_OWNERS = frozenset((
    "unit_test_m1_neutrino_seeded_invariants",
    "unit_test_m1_neutrino_rusanov_flux",
    "unit_test_m1_thcm1_blended_rusanov",
    "unit_test_rusanov_flux",
))
EXPECTED_PRODUCER_STATE = "retained_historical"
EXPECTED_ADMISSION_STATUS = "not_current_verification_admission"
ADMITTED_PRODUCER_STATE = "current_admitted"
ADMITTED_ADMISSION_STATUS = "current_verification_admitted"
ADMISSION_METADATA_KEYS = (
    "canonical_target",
    "artifact_reference",
    "producer_state",
)
HEX_DIGEST = re.compile(r"^[0-9a-f]{64}$")
COMPACT_IDENTIFIER = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]*$")
RECORD_STRING_FIELDS = frozenset((
    "case_id", "pair_id", "origin", "seed_id", "family", "perturbation",
))
RECORD_INTEGER_FIELDS = frozenset((
    "sensitivity_start", "sensitivity_count", "baseline_available",
    "perturbed_available", "baseline_status", "perturbed_status",
    "baseline_reference_status", "perturbed_reference_status",
))
RECORD_VECTOR_FIELDS = frozenset((
    "baseline_input", "perturbed_input", "baseline_output",
    "perturbed_output", "normalization",
))
RECORD_FIELDS = (RECORD_STRING_FIELDS | RECORD_INTEGER_FIELDS |
                 RECORD_VECTOR_FIELDS)


def _reject_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate manifest key: " + str(key))
        result[key] = value
    return result


def _require_keys(value, expected, context):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError("invalid " + context + " keys")


def _require_string(value, context):
    if not isinstance(value, str) or not value or any(character.isspace() for character in value):
        raise ValueError("invalid " + context)
    return value


def _require_compact_identifier(value, context):
    if (not isinstance(value, str) or
            COMPACT_IDENTIFIER.fullmatch(value) is None):
        raise ValueError("invalid " + context)
    return value


def _require_positive_integer(value, context):
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError("invalid " + context)
    return value


def _load_manifest(root):
    manifest_path = root / MANIFEST_NAME
    if not manifest_path.is_file() or manifest_path.is_symlink():
        raise ValueError("missing package manifest: " + str(manifest_path))
    try:
        manifest = json.loads(
            manifest_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_keys,
        )
    except (OSError, UnicodeError, ValueError) as error:
        raise ValueError("malformed package manifest: " + str(error))

    _require_keys(
        manifest,
        ("schema", "version", "digest", "provenance", "payloads"),
        "manifest",
    )
    if manifest["schema"] != SCHEMA or manifest["version"] != VERSION:
        raise ValueError("unsupported package manifest schema or version")

    digest = manifest["digest"]
    _require_keys(digest, ("algorithm", "purpose"), "digest metadata")
    if digest["algorithm"] != "sha256" or digest["purpose"] != "fixture_payload_integrity_only":
        raise ValueError("invalid digest metadata")

    provenance = manifest["provenance"]
    _require_keys(provenance, ("producer_state", "admission_status"), "provenance")
    if (provenance["producer_state"] == EXPECTED_PRODUCER_STATE and
            provenance["admission_status"] == EXPECTED_ADMISSION_STATUS):
        package_is_admitted = False
    elif (provenance["producer_state"] == ADMITTED_PRODUCER_STATE and
          provenance["admission_status"] == ADMITTED_ADMISSION_STATUS):
        package_is_admitted = True
    else:
        raise ValueError("invalid provenance or admission status")

    payloads = manifest["payloads"]
    if not isinstance(payloads, list) or not payloads:
        raise ValueError("manifest payload list is empty or invalid")
    paths = set()
    owners = set()
    for index, payload in enumerate(payloads):
        context = "payload entry " + str(index)
        if package_is_admitted and (
                not isinstance(payload, dict) or "admission" not in payload):
            raise ValueError("missing admission metadata in " + context)
        if (not package_is_admitted and isinstance(payload, dict) and
                "admission" in payload):
            raise ValueError("admission metadata is not allowed for historical package")
        _require_keys(
            payload,
            ("path", "owner", "operation", "policy", "record_count", "sha256") +
            (("admission",) if package_is_admitted else ()),
            context,
        )
        if package_is_admitted:
            admission = payload["admission"]
            _require_keys(admission, ADMISSION_METADATA_KEYS, context + " admission metadata")
            _require_compact_identifier(
                admission["canonical_target"],
                context + " admission metadata canonical target",
            )
            _require_compact_identifier(
                admission["artifact_reference"],
                context + " admission metadata artifact reference",
            )
            if admission["producer_state"] != ADMITTED_PRODUCER_STATE:
                raise ValueError("invalid " + context + " admission metadata producer state")
        path = _require_string(payload["path"], context + " path")
        path_object = Path(path)
        if (path_object.is_absolute() or path_object.name != path or
                ".." in path_object.parts or path_object.suffix not in PAYLOAD_SUFFIXES):
            raise ValueError("invalid " + context + " path")
        if path in paths:
            raise ValueError("duplicate payload path: " + path)
        paths.add(path)
        owner = _require_string(payload["owner"], context + " owner")
        if owner not in STORED_REFERENCE_OWNERS:
            raise ValueError("payload has an unknown stored-reference owner: " + owner)
        owners.add(owner)
        _require_string(payload["operation"], context + " operation")
        _require_string(payload["policy"], context + " policy")
        _require_positive_integer(payload["record_count"], context + " record count")
        if (not isinstance(payload["sha256"], str) or
                HEX_DIGEST.fullmatch(payload["sha256"]) is None):
            raise ValueError("invalid " + context + " payload digest")

    if owners != STORED_REFERENCE_OWNERS:
        raise ValueError("manifest does not identify all stored-reference owners")
    return payloads, paths


def _actual_payload_paths(root):
    actual = set()
    for path in root.iterdir():
        if path.suffix in PAYLOAD_SUFFIXES:
            actual.add(path.name)
    return actual


def _reject_external_sidecars(root):
    for path in root.iterdir():
        if (path.name != MANIFEST_NAME and
                (path.name.endswith(".receipt.json") or
                 path.name.endswith(".manifest.json"))):
            raise ValueError("external campaign sidecar is not public package data: " +
                             path.name)


RECORD_FIELD_ORDER = (
    "case_id", "pair_id", "origin", "seed_id", "family", "perturbation",
    "sensitivity_start", "sensitivity_count", "baseline_input",
    "perturbed_input", "baseline_output", "perturbed_output", "normalization",
    "baseline_available", "perturbed_available", "baseline_status",
    "perturbed_status", "baseline_reference_status", "perturbed_reference_status",
)


def _take(tokens, context):
    try:
        return next(tokens)
    except StopIteration as error:
        raise ValueError("truncated payload: " + context) from error


def _parse_unsigned(token, context):
    if re.fullmatch(r"[0-9]+", token) is None:
        raise ValueError("invalid record integer: " + context)
    return int(token)


def _parse_vector(tokens, context):
    count = _parse_unsigned(_take(tokens, context + " count"), context + " count")
    if count == 0:
        raise ValueError("invalid record vector: " + context)
    values = []
    for index in range(count):
        token = _take(tokens, context + " value " + str(index))
        try:
            value = float(token)
        except ValueError as error:
            raise ValueError("invalid record vector: " + context) from error
        if not math.isfinite(value):
            raise ValueError("non-finite record vector: " + context)
        values.append(value)
    return values


def _validate_record(record, context):
    if set(record) != RECORD_FIELDS:
        raise ValueError("invalid record fields: " + context)
    for field in RECORD_STRING_FIELDS:
        if (not isinstance(record[field], str) or not record[field] or
                any(character.isspace() for character in record[field])):
            raise ValueError("invalid record string: " + context)
    if record["case_id"] == record["pair_id"]:
        raise ValueError("record case and pair IDs collapse: " + context)
    for field in RECORD_INTEGER_FIELDS:
        if not isinstance(record[field], int) or record[field] < 0:
            raise ValueError("invalid record integer: " + context)
    if (record["baseline_available"] != 1 or
            record["perturbed_available"] != 1 or
            record["baseline_status"] != 0 or
            record["perturbed_status"] != 0 or
            record["baseline_reference_status"] != 0 or
            record["perturbed_reference_status"] != 0):
        raise ValueError("record reference status is not successful: " + context)
    baseline_input = record["baseline_input"]
    perturbed_input = record["perturbed_input"]
    baseline_output = record["baseline_output"]
    perturbed_output = record["perturbed_output"]
    normalization = record["normalization"]
    if (len(baseline_input) != len(perturbed_input) or
            len(baseline_output) != len(perturbed_output) or
            len(normalization) != len(baseline_output) or
            any(value <= 0.0 for value in normalization)):
        raise ValueError("record vector shapes are inconsistent: " + context)
    start = record["sensitivity_start"]
    count = record["sensitivity_count"]
    if count == 0:
        if (baseline_input != perturbed_input or
                baseline_output != perturbed_output):
            raise ValueError("input-invariant record has a response: " + context)
        return
    if start >= len(baseline_input) or count > len(baseline_input) - start:
        raise ValueError("record sensitivity range is invalid: " + context)
    if not any(baseline_input[index] != perturbed_input[index]
               for index in range(start, start + count)):
        raise ValueError("record sensitivity range has no perturbation: " + context)
    if any(baseline_input[index] != perturbed_input[index]
           for index in range(len(baseline_input))
           if index < start or index >= start + count):
        raise ValueError("record changes outside sensitivity range: " + context)


def _inspect_payload(path):
    digest = hashlib.sha256()
    def payload_tokens():
        try:
            with path.open("rb") as payload:
                for line_number, raw_line in enumerate(payload):
                    digest.update(raw_line)
                    try:
                        line = raw_line.decode("utf-8")
                    except UnicodeDecodeError as error:
                        raise ValueError(
                            "payload is not UTF-8 at line " + str(line_number + 1)
                        ) from error
                    yield from line.split()
        except OSError as error:
            raise ValueError("cannot read payload: " + str(path)) from error

    try:
        tokens = iter(payload_tokens())
        if _take(tokens, "magic") != "M1_THCM1_FIXTURE":
            raise ValueError("invalid payload header")
        if _parse_unsigned(_take(tokens, "version"), "version") != 1:
            raise ValueError("invalid payload header version")
        if _take(tokens, "operation key") != "operation":
            raise ValueError("invalid payload operation header")
        operation = _take(tokens, "operation")
        if not operation or any(character.isspace() for character in operation):
            raise ValueError("invalid payload operation header")
        if _take(tokens, "policy key") != "policy":
            raise ValueError("invalid payload policy header")
        policy = _take(tokens, "policy")
        if not policy or any(character.isspace() for character in policy):
            raise ValueError("invalid payload policy header")
        if _take(tokens, "record-count key") != "record_count":
            raise ValueError("invalid payload record-count header")
        header_record_count = _parse_unsigned(
            _take(tokens, "record count"), "record count"
        )
        if header_record_count == 0:
            raise ValueError("invalid payload header values")

        case_ids = set()
        pair_ids = set()
        for index in range(header_record_count):
            if _take(tokens, "record begin") != "record_begin":
                raise ValueError("invalid record begin at index " + str(index))
            record = {}
            for field in RECORD_FIELD_ORDER:
                key = _take(tokens, "record field " + field)
                if key != field or key in record:
                    raise ValueError("invalid record field at index " + str(index))
                if field in RECORD_STRING_FIELDS:
                    value = _take(tokens, "record string " + field)
                    if not value:
                        raise ValueError("invalid record string at index " + str(index))
                    record[field] = value
                elif field in RECORD_VECTOR_FIELDS:
                    record[field] = _parse_vector(
                        tokens, str(path) + " record " + str(index) + " " + field
                    )
                else:
                    record[field] = _parse_unsigned(
                        _take(tokens, "record integer " + field),
                        str(path) + " record " + str(index) + " " + field,
                    )
            if _take(tokens, "record end") != "record_end":
                raise ValueError("invalid record termination at index " + str(index))
            _validate_record(record, str(path) + " record " + str(index))
            if (record["case_id"] in case_ids or
                    record["pair_id"] in pair_ids):
                raise ValueError("duplicate record ID at index " + str(index))
            case_ids.add(record["case_id"])
            pair_ids.add(record["pair_id"])

        if _take(tokens, "end marker") != "end":
            raise ValueError("payload is missing end marker")
        try:
            trailing = next(tokens)
        except StopIteration:
            trailing = None
        if trailing is not None:
            raise ValueError("payload has trailing content after end")
    except ValueError:
        raise
    return operation, policy, header_record_count, digest.hexdigest()


def audit(root):
    if not root.is_dir() or root.is_symlink():
        raise ValueError("fixture package root is not a directory")
    _reject_external_sidecars(root)
    payloads, manifest_paths = _load_manifest(root)
    actual_paths = _actual_payload_paths(root)
    missing = sorted(manifest_paths - actual_paths)
    extra = sorted(actual_paths - manifest_paths)
    if missing:
        raise ValueError("missing payload: " + ", ".join(missing))
    if extra:
        raise ValueError("extra payload: " + ", ".join(extra))

    total_records = 0
    for payload in payloads:
        path = root / payload["path"]
        if path.is_symlink() or not path.is_file():
            raise ValueError("payload is not a regular file: " + payload["path"])
        operation, policy, record_count, digest = _inspect_payload(path)
        if (operation != payload["operation"] or policy != payload["policy"] or
                record_count != payload["record_count"]):
            raise ValueError("payload header metadata mismatch: " + payload["path"])
        if digest != payload["sha256"]:
            raise ValueError("payload digest mismatch: " + payload["path"])
        total_records += record_count
    return len(payloads), total_records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    try:
        payload_count, record_count = audit(args.root.resolve())
    except (OSError, ValueError) as error:
        parser.exit(1, "fixture package audit failed: " + str(error) + "\n")
    print("offline M1 fixture package audit passed: " +
          str(payload_count) + " payloads, " + str(record_count) + " records")


if __name__ == "__main__":
    main()
