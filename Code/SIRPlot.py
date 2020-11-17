# defines plot functions

# import libraries
import matplotlib.pyplot as plt
import seaborn as sns
import os

from Code.SIRParams import df

# Force the correct directory
if os.getcwd().split("/")[-1] == "Code":
    os.chdir("..")
curr_dir = os.getcwd()

# If output directory does not already exist, create one
if not os.path.isdir("Output"):
    os.mkdir("Output")
output_dir = curr_dir + "/Output/"

# Plot the data on three separate curves for S(t), I(t) and R(t)

def plot_sir(SIR, params, country):
    SIR_long = SIR.melt(
        id_vars=['Days'], value_vars=['S', 'I', 'R'], value_name='Number of People', var_name='Status')
    fig, ax = plt.subplots(ncols=1, nrows=1)
    sns.lineplot(data=SIR_long, x='Days', y='Number of People', hue='Status', ax=ax)
    ax.set_title('SIR Model ' + country)
    ax.annotate('γ: ' + str(round(params['gamma'], 2)), xy=(1, 0), xycoords='axes fraction', fontsize=12,
                xytext=(-5, 5), textcoords='offset points',
                ha='right', va='bottom')
    fig.set_size_inches(18.5, 10.5)
    plt.savefig(output_dir + "SIR_plot_" + country + ".png")
    plt.close()


def plot_i(SIR, params, country):
    SIR_long = SIR.melt(
        id_vars=['Days'], value_vars=['I'], value_name='Number of People', var_name='Status')
    fig, ax = plt.subplots(ncols=1, nrows=1)
    sns.lineplot(data=SIR_long, x='Days', y='Number of People', hue='Status', ax=ax)
    ax.set_title('Estimated Infected based on SIR ' + country)
    ax.annotate('γ: ' + str(round(params['gamma'], 2)), xy=(1, 0), xycoords='axes fraction', fontsize=12,
                xytext=(-5, 5), textcoords='offset points',
                ha='right', va='bottom')
    fig.set_size_inches(18.5, 10.5)
    plt.savefig(output_dir + "I_plot_" + country + ".png")
    plt.close()


def plot_multiple_sir(SIR, params, country):
    SIR_long = SIR.melt(
        id_vars=['Days', 'R0'], value_vars=['S', 'I', 'R'], value_name='Number of People', var_name='Status')
    fig, ax = plt.subplots(ncols=1, nrows=1)
    sns.lineplot(data=SIR_long, x='Days', y='Number of People', hue='Status', style= 'R0', ax=ax)
    ax.set_title('SIR Model ' + country)
    ax.annotate('γ: ' + str(round(params['gamma'], 2)), xy=(1, 0), xycoords='axes fraction', fontsize=12,
                xytext=(-5, 5), textcoords='offset points',
                ha='right', va='bottom')
    plt.legend(loc='upper right')
    fig.set_size_inches(18.5, 10.5)
    plt.savefig(output_dir + "SIR_plot_multiple_R0_" + country + ".png")
    plt.close()


