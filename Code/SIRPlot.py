# defines plot functions

# import libraries
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Force the correct directory
if os.getcwd().split("/")[-1] == "Code":
    os.chdir("..")
curr_dir = os.getcwd()

# If output directory does not already exist, create one
if not os.path.isdir("Output"):
    os.mkdir("Output")
output_dir = curr_dir + "/Output/"

# Plot the data on three separate curves for S(t), I(t) and R(t)

def plot_sir(SIR, country):
    SIR_long = SIR.melt(
        id_vars=['Days'], value_vars=['S', 'I', 'R'], value_name='Number of People', var_name='Status')
    fig, ax = plt.subplots(ncols=1, nrows=1)
    sns.lineplot(data=SIR_long, x='Days', y='Number of People', hue='Status', ax=ax)
    ax.set_title('SIR Model ' + country)
    plt.savefig(output_dir + "SIR_plot_" + country + ".png")
    plt.close()


