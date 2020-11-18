################################################################################
# Make Mobility Videos
################################################################################

data_dict['mobility']['Date'] = pd.to_datetime(data_dict['mobility']['Date'],
                                               format='%Y-%m-%d')
mob_week = data_dict['mobility'].groupby(
    [data_dict['mobility'].index, 'StateName', pd.Grouper(key='Date', freq='W-MON')]).agg(np.nanmean).reset_index(level=[1, 2])
mob_week = mob_week[~mob_week['StateName'].isin(['AK', 'HI', 'PR', 'VI'])]
mob_week = get_gp_data(mob_week, data_dict['geo_map'])

dates = sorted(mob_week['Date'].unique())

plot_params['col'] = 'Residential'
plot_params['title'] = 'Google Mobility %Change from Baseline - Parks'
plot_params['file'] = output_dir + '/Parks.png'
plot_params['date_bar'] = False
make_choropleth(mob_week[mob_week['Date'] == dates[0]], plot_params)

if not os.path.isdir(col + '_video'):
    os.mkdir(col + '_video')

dates = sorted(mob_week['Date'].unique())

vmin = mob_week[col].min()
vmax = mob_week[col].max()
for i, date in enumerate(dates):
    fig, ax = plt.subplots(1, figsize=(50, 20))
    sub_df = plot_data[plot_data['Date'] == date]
    if sub_df[col].isnull().sum() != 0:
        sub_df.plot(column=col, ax=ax, cmap='viridis',
                    missing_kwds=dict(color='grey', label='No Data'))
    else:
        sub_df.plot(column=col, ax=ax, cmap='viridis')
    ax.axis('off')
    ax.set_title(col + ' Mobility',
                 fontdict={'fontsize': '75', 'fontweight': '15'})
    sm = plt.cm.ScalarMappable(
        cmap='viridis',
        norm=plt.Normalize(vmin=vmin,
                           vmax=vmax))

    cbar = fig.colorbar(sm)
    cbar.ax.tick_params(labelsize=25)

    # Positions for the date
    date_x = -75
    date_y = 50
    ax.text(date_x, date_y,
            np.datetime_as_string(date, unit='D'),
            color='black',
            fontsize=50, fontweight=15)
    fig.savefig(col + f"_video/frame_{i:03d}.png",
                dpi=100, bbox_inches='tight')
    plt.close()

os.chdir(col + "_video")
os.system(
    "ffmpeg -framerate 5 -i frame_%3d.png -c:v h264 -r 30 -s 1920x1080 ./" +
    col + "_mobility.mp4")
os.chdir("..")
