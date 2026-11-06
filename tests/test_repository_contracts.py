from __future__ import annotations

import importlib.util
import tempfile
import unittest
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
P5 = ROOT / "Explainable-AI-Cancer-Drug-Response"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Project5MetricTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rna_models = load_module("project5_rna_models", P5 / "03_rna_models.py")

    def test_mae_and_r2_are_not_swapped(self):
        result = self.rna_models.metrics([0, 1, 2], [0, 1, 1])
        self.assertAlmostEqual(result["MAE"], 1 / 3)
        self.assertAlmostEqual(result["R2"], 0.5)

    def test_established_phase2_output_schema_and_values(self):
        selected = pd.DataFrame({"drug_id": ["A"], "model": ["ElasticNet"], "MAE": [0.08]})
        baseline = pd.DataFrame({"MAE": [0.10]}, index=["A"])
        result = self.rna_models.format_selected_results(selected, baseline)
        self.assertEqual(result.loc[0, "selected_model"], "ElasticNet")
        self.assertAlmostEqual(result.loc[0, "baseline_MAE"], 0.10)
        self.assertAlmostEqual(result.loc[0, "MAE_improvement"], 0.02)
        self.assertTrue(result.loc[0, "improved"])


class Project5RNAReaderTests(unittest.TestCase):
    def test_full_rna_reader_uses_one_identifier_column(self):
        module = load_module("project5_nested_validation", P5 / "06_nested_validation.py")
        model_ids = [f"SIDM{i:04d}" for i in range(300)]
        columns = ["identifier", *model_ids]
        rows = [
            ["model_name", *[f"model_{i}" for i in range(300)]],
            ["data_source", *(["Sanger"] * 300)],
            ["gene_symbol", *([""] * 300)],
            ["GENE1", *np.linspace(1, 2, 300)],
            ["GENE2", *np.linspace(2, 3, 300)],
        ]
        with tempfile.TemporaryDirectory() as directory:
            data_dir = Path(directory)
            csv_path = data_dir / "rna.csv"
            pd.DataFrame(rows, columns=columns).to_csv(csv_path, index=False)
            archive = data_dir / "rnaseq_merged_fixture.zip"
            with zipfile.ZipFile(archive, "w") as handle:
                handle.write(csv_path, arcname="fixture_rsem_tpm_test.csv")
            module.DATA = data_dir
            matrix = module.read_all_rna(set(model_ids))
        self.assertEqual(matrix.shape, (300, 2))
        self.assertEqual(matrix.columns.tolist(), ["GENE1", "GENE2"])


class PipelineContractTests(unittest.TestCase):
    def test_project5_runner_dependency_order_and_files(self):
        runner = load_module("project5_runner", P5 / "run_pipeline.py")
        self.assertEqual(runner.ORDER[:3], ("01_prepare_data.py", "03_rna_models.py", "02_response_landscape.py"))
        self.assertTrue(all(path.is_file() for path in runner.SCRIPTS))

    def test_dynamic_helper_targets_exist(self):
        for relative in ("04_multiomics_explainability.py", "06_nested_validation.py"):
            self.assertTrue((P5 / relative).is_file(), relative)

    def test_all_five_projects_have_ordered_runners(self):
        runners = (
            ROOT / "Single-Cell-RNAseq-Lung-Disease" / "run_pipeline.R",
            ROOT / "Spatial-Transcriptomics-Breast-Cancer" / "run_pipeline.py",
            ROOT / "Pan-Cancer-Transcriptomic-Survival" / "run_pipeline.R",
            ROOT / "Shotgun-Metagenomics-Colorectal-Cancer" / "run_pipeline.R",
            P5 / "run_pipeline.py",
        )
        self.assertTrue(all(path.is_file() for path in runners))


if __name__ == "__main__":
    unittest.main()
