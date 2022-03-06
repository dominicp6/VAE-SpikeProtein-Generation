
def get_delta_delta_g_from_file(filename):

    total_delta_delta_g = 0

    with open(filename, 'w') as file:
        for line_number, line in enumerate(file):
            if line_number == 0:
                continue
            else:
                line_segments = line.split(' ')
                total_delta_delta_g += float(line_segments[-1])

    return total_delta_delta_g


if __name__ == "__main__":

    for seq_number in range(141):
        pass
