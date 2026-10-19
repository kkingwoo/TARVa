"""
The p-value convention every analysis step shares, in one place.

It lives here rather than being reimplemented in each of the six steps that
report a p-value, so the same decision cannot be made differently in two of
them.

The cascade, applied the same way at every step:

  1. BH-adjust. If at least MIN_SIGNIFICANT_FRACTION of the step's tests
     clear alpha on the ADJUSTED values, those are the result.
  2. Otherwise say so, and if that fraction clear alpha on the NOMINAL (raw)
     values, those are the result instead.
  3. Otherwise neither basis qualifies. Nothing is significant, and the
     message says what a user can do about it: more individuals per
     condition, a different dataset, or a lower --min-depth.

Outcome 3 is the one that matters. A step that fell back to nominal p-values
and kept whatever cleared -- three genes out of two thousand, say -- would
return a gene set that is noise wearing the name of a result. The fraction
requirement stops that, and when it is not met the output is nothing plus an
explanation.

`purpose` does not change the cascade, only the wording and what an empty
result means downstream:

  GATING     -- the step sorts genes or sites into categories a later step
                depends on, so an empty result selects no genes.
  REPORTING  -- the p-values are themselves the result, so an empty result
                means nothing is presented as significant.

Input:  nothing -- this module is imported by the analysis steps.
Output: the chosen p-value basis and the message explaining it.
"""

import numpy as np

from tarva.dbpipe import schema

GATING, REPORTING = "gating", "reporting"

ADJUSTED, RAW = "adj_pvalue", "p_value"

# A basis is adopted only when at least this fraction of the step's tests
# clear alpha on it. Adjusted is checked first, nominal only if adjusted
# fails, and neither qualifying means nothing is significant.
#
# Not a majority: 30% is enough to treat a basis as supported. The bar is
# deliberately below half because a real filtering step can leave most genes
# non-significant and still have selected something meaningful -- what it
# cannot do is select a handful of borderline tests and call them a gene set,
# which is what the third outcome exists to prevent.
MIN_SIGNIFICANT_FRACTION = 0.30


P_DECIMALS = 3


def format_p(value, places=P_DECIMALS):
    """
    A p-value at three decimal places -- except where that would print 0.000,
    which is shown in exponent form instead.

    Three decimals is the reporting convention, and it is also roughly where a
    small p-value stops being readable. Flattening 2.448e-160 to 0.000 would
    throw away the only interesting thing about it, so anything below half of
    the last retained place keeps three significant digits and an exponent.

    DISPLAY ONLY. Nothing stored is rounded: gene_test, gene_calltype,
    calltype_test, run_meta and every CSV carry full float precision, and
    to_sql / to_csv are called without float_format so the R scripts threshold
    on untruncated values.

    The rounding is ONE step, taken at the fourth decimal against the whole
    remaining value -- not the fifth into the fourth and then the fourth into
    the third. 0.09446 prints 0.094, never 0.095. Carries propagate as far as
    they need to, so 0.0996 prints 0.100. (Python breaks an exact .0005 tie to
    the even digit; a p-value is never exactly that in binary, so this is
    theoretical.)

    NaN prints as "na": not testable is a different statement from 1.000.
    """
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "na"
    if np.isnan(value):
        return "na"
    if value and abs(value) < 0.5 * 10 ** -places:
        return f"{value:.{places}e}"
    return f"{value:.{places}f}"


def benjamini_hochberg(p_values):
    """
    BH step-up adjusted p-values, input order preserved, NaN passed through.

    NaN means "not testable" -- both arms all-zero, say -- which is a
    different statement from a p-value of 1, so those rows are excluded from
    the ranking rather than treated as non-significant.
    """
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan)
    testable = ~np.isnan(values)
    if not testable.any():
        return adjusted

    present = values[testable]
    order = np.argsort(present)
    ranked = present[order]
    n = ranked.size
    scaled = ranked * n / np.arange(1, n + 1)
    # Step up from the largest p so the sequence is monotone.
    monotone = np.minimum.accumulate(scaled[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(monotone, 1.0)
    adjusted[testable] = out
    return adjusted


def resolve_gene_set(selected, available, label, source):
    """
    The genes a step operates on: the filtered set when it has members,
    otherwise every gene available to it.

    A filter that selects nothing would stop every later step, so the fallback
    keeps the pipeline runnable end to end. It warns, because a result computed
    over all genes is NOT a result computed over a filtered set and must not be
    reported as one -- run_meta records which it was.

    This is a supported outcome, not a diagnostic. A dataset that legitimately
    yields no global gene set still gets every downstream result, over the
    genes that passed the entry criterion, with the message saying which set
    the numbers describe.

    Note that callers filter on the returned set either way -- in the fallback
    it simply contains everything available -- so neither branch takes a
    different code path downstream.

    Returns (genes, filtered); `filtered` is False when the fallback fired.
    """
    if selected:
        return set(selected), True
    schema.warn(f"{label}: {source} selected no genes, so this step runs on ALL "
                f"{len(available):,} gene(s) that passed the entry criterion instead. "
                f"The result is real but it describes the UNFILTERED gene set, not a "
                f"{source} subset, and must be reported that way. run_meta records "
                f"which set was used.")
    return set(available), False


def choose_basis(tests, alpha, purpose, label, downstream=None):
    """
    Which p-value column this step's significance is read from, and how many
    clear it.

    `tests` is a frame already carrying both ADJUSTED and RAW columns. Returns
    (column, n_significant), where the column is None when neither basis
    reaches alpha in enough of the tests. Counts are reported either way, so
    an empty result is never unexplained.

    `downstream` names what cannot run if nothing is significant -- it only
    affects the message.
    """
    if purpose not in (GATING, REPORTING):
        raise ValueError(f"purpose must be {GATING!r} or {REPORTING!r}, not {purpose!r}")

    total = len(tests)
    n_adjusted = int((tests[ADJUSTED] < alpha).sum())
    n_raw = int((tests[RAW] < alpha).sum())
    adjusted_share = n_adjusted / total if total else 0.0
    raw_share = n_raw / total if total else 0.0

    schema.report(f"{label}: {n_adjusted:,} of {total:,} at adjusted p < {alpha} "
                  f"({adjusted_share:.1%}), {n_raw:,} at nominal p < {alpha} "
                  f"({raw_share:.1%})")

    if adjusted_share >= MIN_SIGNIFICANT_FRACTION:
        schema.report(f"{label}: adjusted p-values reach {alpha} in at least "
                      f"{MIN_SIGNIFICANT_FRACTION:.0%} of tests, so this step reports "
                      f"{n_adjusted:,} ADJUSTED result(s)")
        return ADJUSTED, n_adjusted

    schema.report(f"{label}: under {MIN_SIGNIFICANT_FRACTION:.0%} of tests survive "
                  f"correction, so the adjusted p-values do not support a result here; "
                  f"checking nominal")

    if raw_share >= MIN_SIGNIFICANT_FRACTION:
        schema.report(f"{label}: nominal p-values reach {alpha} in at least "
                      f"{MIN_SIGNIFICANT_FRACTION:.0%} of tests, so this step uses "
                      f"NOMINAL p < {alpha}, keeping {n_raw:,}")
        return RAW, n_raw

    # Neither basis qualifies. Keeping whichever genes happened to clear a
    # nominal threshold would name a handful of borderline tests a "global
    # gene set", so the step selects nothing and says why.
    message = (f"{label}: neither basis reaches {alpha} in {MIN_SIGNIFICANT_FRACTION:.0%} "
               f"of the {total:,} test(s) -- {n_adjusted:,} adjusted, {n_raw:,} nominal "
               f"-- so nothing is significant and this step selects nothing.")
    if downstream:
        message += f" {downstream} has nothing to work on."
    message += (f" The test statistics, effect sizes and p-values all remain "
                f"reportable. If a result is expected here, the usual causes are too "
                f"few individuals per condition, a dataset without the signal, or a "
                f"--min-depth strict enough to have excluded the supporting calls.")
    schema.warn(message)
    return None, 0
