import os
import sys
import pandas as pd

TEST_DIR = os.path.dirname(__file__)
RESULTS_DIR = os.path.join(TEST_DIR, "results")
if not os.access(RESULTS_DIR, os.W_OK):
    # repo dir is currently read-only during github-actions testing
    RESULTS_DIR = '/tmp/unpast/results'
REFERENCE_OUTPUT_DIR = os.path.join(TEST_DIR, "test_reference_output")

sys.path.append(os.path.join(TEST_DIR, ".."))
from run_unpast import run

### Helper functions ###

def run_unpast_on_file(filename, basename):
    run(
        os.path.join(TEST_DIR, filename),
        out_dir=RESULTS_DIR,
        basename=basename,
    )
    return parse_answer(RESULTS_DIR, basename)


def parse_answer(answer_dir, basename):
    files = os.listdir(answer_dir)
    answer_files = [f for f in files if f.startswith(basename) and f.endswith("biclusters.tsv")]
    assert len(answer_files) == 1, f"There are {len(answer_files)} files instead of 1"
    return pd.read_csv(os.path.join(answer_dir, answer_files[0]), sep="\t", comment='#')

### Tests ###

def test_smoke():
    """Smoke test - check that the program runs on some input without failure."""
    run_unpast_on_file(
        filename = "test_input/synthetic_clear_biclusters.tsv",
        basename = "test_smoke",
    )


def test_clear_biclusters():
    """Check that clear biclusters are found."""
    res = run_unpast_on_file(
        filename = "test_input/synthetic_clear_biclusters.tsv",
        basename = "test_clear_biclusters",
    )
    raise NotImplementedError("TODO: Implement this test")


def test_reproducible():
    """Check that the same data is found on a complicated input with no clear answer."""
    res = run_unpast_on_file(
        filename = "test_input/synthetic_noise.tsv",
        basename = "test_reproducible",
    )
    reference = parse_answer(
        answer_dir = REFERENCE_OUTPUT_DIR,
        basename = "test_reproducible",
    )
    assert res.equals(reference), "The results are not reproducible"


if __name__ == "__main__":
    # run all the tests in this file
    test_smoke()
    test_reproducible()

    # deselected test
    # test_clear_biclusters()

    # TODO: use pytest instead
    # pytest.main([f"{os.path.abspath(__file__)}"])
