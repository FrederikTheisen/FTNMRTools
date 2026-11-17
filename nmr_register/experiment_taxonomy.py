"""Shared experiment names, broad types, and conservative pulse-name recognition.

Pulse names are evidence for classification only; callers retain the raw name.
Unknown explicit experiment names survive normalization and belong to Other.
"""
import re
import inspect


def _key(value):
    return re.sub(r"[\s_-]+", "", value).casefold()


_GROUPS = {
    "Relaxation": ("T1", "T2", "T1rho", "HetNOE", "etaXY", "Water T1 (T1sat)", "CPMG"),
    "HSQC": ("HSQC", "13C HSQC"),
    "Assignment": ("HNcaCO", "HNCA", "HNCO", "HNCACB", "HNcoCACB", "HNcoCA", "HNCANH", "iHNCACB", "HNN"),
    "Calibration": ("Calibration", "15N calibration", "13C calibration", "zg", "zg30"),
    "1D": ("1D",),
    "Other": ("HMBC", "COSY", "NOESY"),
}
_ALIASES = {_key(name): name for names in _GROUPS.values() for name in names}
_ALIASES.update({
    "hn(ca)co": "HNcaCO", "hn(co)ca": "HNcoCA", "hn(co)cacb": "HNcoCACB",
    "r1": "T1", "15nt1": "T1", "r2": "T2", "r1rho": "T1rho", "t1ρ": "T1rho",
    "noe": "HetNOE", "heteronuclearnoe": "HetNOE", "ηxy": "etaXY",
    "15nhsqc": "HSQC", "nhsqc": "HSQC", "watert1": "Water T1 (T1sat)",
    "watert1sat": "Water T1 (T1sat)", "t1sat": "Water T1 (T1sat)",
    "nicot1sat": "Water T1 (T1sat)", "15ncalib": "15N calibration", "13ccalib": "13C calibration",
})


def canonical_experiment(value):
    """Return a canonical precise label, retaining unfamiliar explicit labels."""
    if not isinstance(value, str) or not value.strip() or value.strip().casefold() == "unknown":
        return None
    return _ALIASES.get(_key(value), value.strip())


def category(value, pulse=""):
    """Return the broad type without discarding an unfamiliar experiment."""
    name = canonical_experiment(value)
    if not name:
        name = classify_pulse(pulse)
    return next((group for group, names in _GROUPS.items() if name in names), "Other")


def is_calibration(value, pulse=""):
    return category(value, pulse) == "Calibration"


def classify_pulse(pulse, indirect_nucleus=""):
    """Recognize reviewed local naming families; unfamiliar names stay unknown."""
    lower = str(pulse or "").strip().casefold()
    calibration = {"zg": "zg", "zg30": "zg30", "m.zg": "zg",
        "15n_calib": "15N calibration", "13c_calib": "13C calibration",
        "calib15n_fast": "15N calibration", "calib13c_fast": "13C calibration",
        "m.calib15n.wg": "15N calibration"}
    if lower in calibration:
        return calibration[lower]
    # Restrict assignment matching to complete known local families. Prefixes or
    # arbitrary substring matches alone cannot establish a triple-resonance type.
    match = re.fullmatch(r"(?:b_|bt_)(ihncacb|hncocacb|hncoca|hncaco|hncacb|hncanh|hnca|hnco)\+?", lower)
    if match:
        return {"hncaco": "HNcaCO", "hncocacb": "HNcoCACB", "hncoca": "HNcoCA", "ihncacb": "iHNCACB"}.get(match[1], match[1].upper())
    if lower == "ibs_1d_sculpting" or lower in ("presat", "zgpr", "zggpw5", "m.zggpw5"):
        return "1D"
    for token, name in (("hetnoe", "HetNOE"), ("hsqcnoe", "HetNOE"), ("nnh_tr", "etaXY"),
                        ("t1sat", "Water T1 (T1sat)"), ("cpmg", "CPMG"), ("rexet", "CPMG"),
                        ("t1rho", "T1rho"), ("hsqctret", "T1rho"), ("t2", "T2"), ("t1", "T1")):
        if token in lower:
            return name
    if "hsqc" in lower or lower in ("btrosy", "n_sofast_trosy", "n_sofast"):
        return "13C HSQC" if str(indirect_nucleus).strip("<> ").upper() == "13C" else "HSQC"
    for token in ("hmbc", "cosy", "noesy"):
        if token in lower:
            return token.upper()
    return None


EXPERIMENT_CHOICES_BY_CATEGORY = {group: list(names) for group, names in _GROUPS.items()}
EXPERIMENT_CATEGORIES = tuple(_GROUPS)


def choices_by_category():
    return {group: list(names) for group, names in _GROUPS.items()}


experiment_choices_by_category = choices_by_category


def category_formula(value_expr='$ExperimentType', pulse_expr='$PulseProgram'):
    """Build a self-contained Grist formula from this module's canonical aliases."""
    helpers = '\n\n'.join(inspect.getsource(helper) for helper in
                            (_key, canonical_experiment, classify_pulse, category))
    return ("import re\n_GROUPS = " + repr(_GROUPS) + "\n_ALIASES = " + repr(_ALIASES)
            + "\n" + helpers + "\nreturn category(" + value_expr + ", " + pulse_expr + ")")
