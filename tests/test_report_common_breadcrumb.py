from report_common import breadcrumb_nav_html


def test_breadcrumb_runs_broadest_to_most_specific():
    html = breadcrumb_nav_html()
    assert (html.index("All studies") < html.index("Group")
            < html.index("Study</a>") < html.index("Run</a>"))


def test_breadcrumb_paths_match_per_run_layout():
    # Pages sit in docs/<domain>/<study>/<run>/.
    html = breadcrumb_nav_html()
    assert '<a href="../../../index.html">&#127968; All studies</a>' in html
    assert '<a href="../../index.html">&#128194; Group</a>' in html
    assert '<a href="../report.html">&#128193; Study</a>' in html
    assert '<a href="report.html">&#128196; Run</a>' in html


def test_breadcrumb_without_run_link():
    html = breadcrumb_nav_html(run=False)
    assert 'href="report.html"' not in html
    assert "Study</a>" in html
