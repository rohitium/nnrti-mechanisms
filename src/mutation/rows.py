from __future__ import annotations


def build_rows(results, wt_metrics):
    rows = []
    for result in results:
        base = result["base"]
        mut_metrics = result["mut_metrics"]
        for state, metrics in [("WT", wt_metrics), ("MUT", mut_metrics)]:
            rows.extend(
                [
                    {
                        **base,
                        "state": state,
                        "metric": "binding_delta_g_kj_mol",
                        "value": metrics["binding_delta_g_kj_mol"],
                    },
                    {
                        **base,
                        "state": state,
                        "metric": "binding_proxy_kj_mol",
                        "value": metrics["binding_proxy_kj_mol"],
                    },
                    {
                        **base,
                        "state": state,
                        "metric": "contact_count",
                        "value": metrics["contact_count"],
                    },
                    {
                        **base,
                        "state": state,
                        "metric": "hbond_count",
                        "value": metrics["hbond_count"]
                        if metrics["hbond_count"] is not None
                        else float("nan"),
                    },
                    {
                        **base,
                        "state": state,
                        "metric": "pocket_volume_proxy",
                        "value": metrics["pocket_volume_proxy"],
                    },
                ]
            )
    return rows
