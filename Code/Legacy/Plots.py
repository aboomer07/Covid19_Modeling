import matplotlib
from matplotlib import cm
import plotly
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import plotly.express as px
import datetime


def make_choropleth(data, params):
    data = data.copy(deep=True)
    params = params.copy()

    size = params['size']
    col = params['col']
    cmap = params['cmap']
    legend = params['legend']
    title = params['title']
    title_font = params['title_font']
    bar_font = params['bar_font']
    date_bar = params['date_bar']
    file = params['file']

    fig, ax = plt.subplots(1, figsize=size)
    if data[col].isnull().sum() != 0:
        data.plot(column=col, ax=ax, cmap=cmap, legend=legend,
                  missing_kwds=dict(color='lightgrey'))
    else:
        data.plot(column=col, ax=ax, cmap='viridis', legend=legend)

    ax.axis('off')
    ax.set_title(title, fontdict=title_font)
    colorbar = ax.get_figure().get_axes()[1]
    yticks = list(colorbar.get_yticks())

    if date_bar:
        epoch = datetime.datetime(1970, 1, 1)

        result = False
        while not result:
            try:
                colorbar.set_yticks(yticks)
                colorbar.set_yticklabels(
                    [(epoch + datetime.timedelta(seconds=ytick)).
                        strftime('%Y-%m-%d') for ytick in yticks],
                    fontdict=bar_font)
                result = True
            except:
                pass
    else:
        colorbar.set_yticks(yticks)
        colorbar.set_yticklabels(yticks, fontdict=bar_font)

    fig.savefig(file, dpi=100, bbox_inches='tight')
    plt.close()
