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

def read_matrix_from_binary(filename):
    file = open(filename, "rb")
        
    header_format = "qqq"  # Format for long long (n, m, d)
    alignment = 8  # Assume double alignment (sizeof(double) = 8 bytes)
    
    header_format = "qqq"  # Format for long long (n, m, d)
    # Serialize the header
    header_size = struct.calcsize(header_format)
    # Read the first 24 bytes (3 long long values, 8 bytes each)
    data = file.read(24)
    
    if len(data) < 24:
        raise ValueError("File is too small to contain 3 long long numbers.")
    
    # Unpack the data as three long long integers ('q' format specifier)
    n, m, d = struct.unpack(header_format, data)
    alignment = 8  # Assume double alignment (sizeof(double) = 8 bytes)
    
    # Serialize the header
    header_size = struct.calcsize(header_format)
    aligned_header_size = round_up(header_size, alignment)
    padding_length = aligned_header_size - header_size
    file.read(padding_length)

    # Calculate the size of the matrix
    matrix_size = n * n * 8  # Each double is 8 bytes

    # Read the matrix data
    matrix_data = file.read(matrix_size)
    file.close()
    
    if len(matrix_data) < matrix_size:
        raise ValueError("File does not contain enough data for the matrix.")
    
    # Convert the binary data to a numpy array
    X = np.frombuffer(matrix_data, dtype=np.float64).reshape((n, n))
    
    return X
