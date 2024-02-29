import math
def best_layout(num_subplots):
    """
    Determine the best layout for a given number of subplots.

    Parameters:
        - num_subplots: The total number of subplots

    Returns:
        - rows, columns: Best number of rows and columns for the layout
    """

    # If it's a perfect square
    if math.sqrt(num_subplots).is_integer():
        return int(math.sqrt(num_subplots)), int(math.sqrt(num_subplots))

    # If not a perfect square
    for i in range(int(math.sqrt(num_subplots)), 0, -1):
        if num_subplots % i == 0:
            return i, num_subplots // i

    raise ValueError("Unable to determine layout")