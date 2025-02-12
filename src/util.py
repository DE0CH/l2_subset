import struct
import numpy as np

def round_up(location, alignment):
    return (location + alignment - 1) & ~(alignment - 1)

def save_matrix_to_binary(filename, X, n, m):
    
    if m > n:
        raise ValueError("m cannot be greater than n.")
    
    # Use the first batch of X

    header_format = "qqq"  # Format for long long (n, m, d)
    alignment = 8  # Assume double alignment (sizeof(double) = 8 bytes)
    
    # Serialize the header
    header_size = struct.calcsize(header_format)
    aligned_header_size = round_up(header_size, alignment)
    padding_length = aligned_header_size - header_size
    header = struct.pack(header_format, n, m, 0) + b'\x00' * padding_length # d is 0 because we don't really care about
    # Write to file
    with open(filename, "wb") as f:
        f.write(header)
        # we don't write the points array because it's a 0 length array as d = 0
        X.tofile(f)