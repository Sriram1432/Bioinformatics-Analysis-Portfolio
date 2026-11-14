# Explainable AI for cancer drug response

This workflow models drug sensitivity from transcriptomic and multi-omics features, evaluates held-out and nested-validation performance, examines cross-platform transfer and reports SHAP-based model explanations. Predictive associations and feature attributions are not interpreted as causal effects.

Inspect the dependency order without running the analyses:

```text
python run_pipeline.py --dry-run
```

Run the full workflow with `python run_pipeline.py`. The RNA modelling script intentionally runs before the response-landscape script because the latter can consume its held-out predictions. Source datasets and generated model outputs are excluded and are required for numerical reproduction.
