import numpy as np
import matplotlib.pyplot as plt

def length_distribution_single(fig,ax, sequencing_report, sample, font_settings):
    batch = sequencing_report[sequencing_report["Experiment"] == sample]
    length = batch["aaSeqCDR3"].str.len()
    unique_length, counts_length = np.unique(np.array(length)
                                                , return_counts = True)
    ax.bar(unique_length, counts_length)  # Or whatever you want in the subplot
    #ax.set_xticks(range(0, max_length + 1, 1), range(0, max_length + 1, 1))
    ax.title.set_text(sample)
    ax.title.set_size(10)
    ax.set_ylabel("Read Count",
                    **font_settings)  # Y label
    ax.set_xlabel('Read Length',
                **font_settings)  # X label

    original_fontsize = font_settings["fontsize"]
    font_settings["fontsize"] = 22
    plt.title("Length Distribution of " + sample,
              pad=12,
              **font_settings)
    font_settings["fontsize"] = original_fontsize



def length_distribution_multi(fig, sequencing_report, samples,num_cols, font_settings, test_version = False):
    if samples == "all":
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)
    else:
        sequencing_report = sequencing_report[sequencing_report['Experiment'].isin(samples)]
        unique_experiments = sequencing_report["Experiment"].unique()
        unique_experiments = np.sort(unique_experiments)

    # Subplots are organized in a Rows x Cols Grid
    # Tot and Cols are known
    Tot = unique_experiments.shape[0]
    Cols = num_cols
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
       # ax.set_xticks(range(0, max_length + 1, 1), range(0, max_length + 1, 1))
        ax.title.set_text(experiment)
        ax.title.set_size(10)

        adapted_fontsize = 10 - int(Cols) + 2
        font_settings["fontsize"] = adapted_fontsize
        ax.set_ylabel("Read Count",
                      **font_settings)  # Y label
        ax.set_xlabel('Read Length',
                      **font_settings)  # X label

        n += 1
    fig.suptitle("Length Distribution of given Samples")







 #   plt.show()
   # plt.bar(unique_length, counts_length)
   # p_values = counts_length/np.sum(counts_length)
    #normalized on one:
   # plt.bar(unique_length, counts_length/np.max(counts_length))
    #p-values
  #  plt.bar(unique_length, p_values)