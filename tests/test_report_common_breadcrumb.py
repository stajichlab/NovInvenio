from report_common import breadcrumb_nav_html


def test_breadcrumb_runs_broadest_to_most_specific():
    html = breadcrumb_nav_html()
    assert html.index("All studies") < html.index("Group") < html.index("Study</a>")


def test_breadcrumb_without_study_link():
    html = breadcrumb_nav_html(study=False)
    assert "report.html" not in html
    assert html.index("All studies") < html.index("Group")
