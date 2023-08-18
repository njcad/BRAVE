
import numpy as np
from scipy.optimize import linear_sum_assignment
import csv

# csv file with BRAVE and Pilot data
DATA = "/Users/natecadicamo/Desktop/BioX Project/total_demos.csv"

# male controls – note that C008 had no fMRI data, so was removed
controls = {'C002': 69.21, 'C003': 78.05, 'C004': 73.03, 'C005': 64.81, 'C006': 60.1, 'C007': 49.07,
            'C009': 35.22, 'C010': 73.81, 'C011': 51.29}


def build_dicts():
    """
    From BRAVE and Pilot data stored in total_demos.csv, sort subjects into abstainer or relapser dictionary.
    Only make dictionary for males, since we only have 1 female in control group, and have her matches.
    """

    # initialize dictionary of structure {SUBID: age}
    abstainers = {}
    relapsers = {}

    # perform search and build
    with open(DATA, 'r') as f:
        reader = csv.reader(f)
        # skip header row
        next(reader)
        # iterate over all the relevant data
        for row in reader:
            # skip over tagged participants with missing data
            if 'X' in row[0]:
                pass
            # otherwise, they should have data
            else:
                sex = row[3]
                # only consider male participants
                if sex == '2' or sex == 'Male':
                    relapse = row[1]
                    if relapse == '1':
                        # if they relapsed and age is valid, add them to relapsers
                        try:
                            relapsers[row[0]] = float(row[2])
                        except ValueError:
                            # otherwise, skip them
                            pass
                    elif relapse == '0':
                        # if they abstained and age is valid, add them to abstainers
                        try:
                            abstainers[row[0]] = float(row[2])
                        except ValueError:
                            # otherwise, skip them
                            pass

    # return the complete dictionaries
    return abstainers, relapsers


def match_keys_minimize_difference(dict1, dict2):
    """
    Inputs: dict1 and dict2, which will be the male controls and male abstainers / relapsers.
    Outputs: hungarian algorithm derived SUBID matches, and the average difference.
    """
    keys1, values1 = zip(*dict1.items())
    keys2, values2 = zip(*dict2.items())

    # create a cost matrix containing absolute differences between values
    cost_matrix = np.abs(np.subtract.outer(values1, values2))

    # find the optimal pairing using the Hungarian algorithm
    row_indices, col_indices = linear_sum_assignment(cost_matrix)

    # create a list of matched key pairs based on the optimal pairing
    matched_keys = [(keys1[i], keys2[j]) for i, j in zip(row_indices, col_indices)]

    # calculate the total difference and the average difference
    total_diff = sum(cost_matrix[row_indices, col_indices])
    average_diff = total_diff / len(matched_keys)

    return matched_keys, average_diff


def main():
    # get abstainers and relapsers dictionaries
    a, r = build_dicts()

    # match controls to abstainers
    controls_abstainers, ca_delta = match_keys_minimize_difference(controls, a)

    # match controls to relapsers
    controls_relapsers, cr_delta = match_keys_minimize_difference(controls, r)

    # display results
    print("Controls-Abstainers Data:")
    print(f'Matches: {controls_abstainers}')
    print(f'Average difference in age: {ca_delta}\n')
    print("Controls-Relapsers Data:")
    print(f'Matches: {controls_relapsers}')
    print(f'Average differnce in age: {cr_delta}')


if __name__ == '__main__':
    main()