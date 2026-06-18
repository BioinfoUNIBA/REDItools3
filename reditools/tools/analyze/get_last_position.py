import os

def get_pos(redi_bytes_line):
    if redi_bytes_line is None:
        return None
    fields = redi_bytes_line.decode('utf-8').split('\t')
    if len(fields) == 14:
        return int(fields[1])
    return None

def get_last_position(filename):
    with open(filename, 'rb') as stream:
        file_size = os.path.getsize(filename)
        stream.seek(max(0, file_size - 1000))

        next_to_last_line = next(stream, None)
        if next_to_last_line is None:
            return 0

        last_line = next(stream, None)
        for line in stream:
            next_to_last_line = last_line
            last_line = line

        last_pos = get_pos(last_line) 
        if last_pos is not None:
            return last_pos
        last_pos = get_pos(next_to_last_line)
        if last_pos is None:
            return 0
        return last_pos
