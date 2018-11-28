import nbformat as nbf
nb = nbf.v4.new_notebook()

text = """\
# My first automatic Jupyter Notebook
This is an auto-generated notebook."""

code = """\
%pylab inline
hist(normal(size=2000), bins=50);"""

nb['cells'] = [nbf.v4.new_markdown_cell(text),
               nbf.v4.new_code_cell(code) ]

nbf.write(nb, 'TriPlots.ipynb')
