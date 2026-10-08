from src.pathway import prepare_ranked_gene_list


def test_ranked_gene_list_is_sorted(de_results):

    ranked = prepare_ranked_gene_list(
        de_results,
        gene_column="symbol",
        score_column="statistic"
    )

    assert len(ranked) > 0