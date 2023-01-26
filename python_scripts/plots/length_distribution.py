import numpy as np
import matplotlib.pyplot as plt

def length_distribution(fig, sequencing_report, samples, font_settings):
    if samples == True:
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)
    else:
        sequencing_report = sequencing_report[sequencing_report['Experiment'].isin(samples)]
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)

    # Subplots are organized in a Rows x Cols Grid
    # Tot and Cols are known
    Tot = unique_experiments.shape[0]
    Cols = int(input("How many Rows for the Plots do you want?"))
    # Compute Rows required
    Rows = Tot // Cols
    #     EDIT for correct number of rows:
    #     If one additional row is necessary -> add one:
    if Tot % Cols != 0:
        Rows += 1
    # Create a Position index
    Position = range(1,Tot + 1)
    n = 0
   # fig = plt.figure(1, constrained_layout=True)
    for experiment in unique_experiments:
        batch = sequencing_report[sequencing_report["Experiment"] == experiment]
        length = batch["aaSeqCDR3"].str.len()
        unique_length, counts_length = np.unique(np.array(length)
                                                 , return_counts = True)

            # add every single subplot to the figure with a for loop
        ax = fig.add_subplot(Rows, Cols, Position[n])
        ax.bar(unique_length, counts_length)  # Or whatever you want in the subplot
        ax.title.set_text(experiment)
        ax.title.set_size(10)
        adapted_fontsize = 10 - Cols + 2
        font_settings["fontsize"] = adapted_fontsize
        ax.set_ylabel("Read Count",
                      **font_settings)  # Y label
        ax.set_xlabel('Read Length',
                      **font_settings)  # X label
        n += 1
    fig.suptitle('Distribution of number of sequences with certain length')
    plt.show()





 #   plt.show()
   # plt.bar(unique_length, counts_length)
   # p_values = counts_length/np.sum(counts_length)
    #normalized on one:
   # plt.bar(unique_length, counts_length/np.max(counts_length))
    #p-values
  #  plt.bar(unique_length, p_values)