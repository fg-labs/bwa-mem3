"""Keep replacement MSan checks visible without running light-PR workloads."""

from pathlib import Path
import re
import unittest


class MsanCheckLifecycleTest(unittest.TestCase):
    """Check scheduling and step gates in the MSan workflow block."""

    @classmethod
    def setUpClass(cls) -> None:
        workflow = Path(__file__).resolve().parents[1] / "workflows" / "ci.yml"
        text = workflow.read_text()
        block = re.search(
            r"^  msan-padded-lane-gate:\n(.*?)(?=^  [a-z][a-z0-9-]*:|\Z)",
            text,
            re.MULTILINE | re.DOTALL,
        )
        assert block is not None, "MSan job is missing"
        cls.block = block.group(1)

    def test_check_has_no_job_dependency_or_depth_condition(self) -> None:
        """The queued check must exist before any other runner is available."""
        self.assertNotRegex(self.block, re.compile(r"^    needs:", re.MULTILINE))
        condition = re.search(r"^    if: (.+)$", self.block, re.MULTILINE)
        self.assertIsNotNone(condition)
        assert condition is not None
        self.assertNotIn("depth", condition.group(1))
        self.assertNotIn("needs", condition.group(1))

    def test_depth_is_local_and_only_expensive_steps_are_gated(self) -> None:
        """Light PRs skip work, while full or missing depth still runs it."""
        steps = dict(
            re.findall(
                r"^      - name: ([^\n]+)\n(.*?)(?=^      - name:|\Z)",
                self.block,
                re.MULTILINE | re.DOTALL,
            )
        )
        decision = steps["Decide the CI depth"]
        self.assertIn("id: depth", decision)
        self.assertIn("run: bash .github/scripts/ci-depth.sh", decision)
        self.assertNotIn("if:", decision)
        expensive = (
            "Install dependencies (Linux)",
            "Install clang-19 (Linux)",
            "Set up ccache",
            "Build the padded-lane gate under MemorySanitizer",
            "Run the padded-lane gate",
        )
        for name in expensive:
            with self.subTest(step=name):
                self.assertIn("if: steps.depth.outputs.depth != 'light'", steps[name])
        self.assertIn(
            "if: steps.depth.outputs.depth == 'light'",
            steps["Skip the padded-lane workload on a light PR"],
        )


if __name__ == "__main__":
    unittest.main()
