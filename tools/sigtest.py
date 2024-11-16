import signal
import time
import os

def handle_signal(signum, frame):
    elapsed_time = time.time() - start_time
    print(f"Caught signal {signum} ({signal.Signals(signum).name}) at elapsed time {elapsed_time:.1f} seconds")
    print("This program is designed to ignore all signals. To force kill, use `kill -9 {}`".format(os.getpid()))

# Register all signal handlers
for i in range(1, signal.NSIG):
    try:
        signal.signal(i, handle_signal)
    except (OSError, RuntimeError):
        # Some signals cannot be caught or are invalid
        continue

start_time = time.time()
while True:
    elapsed_time = time.time() - start_time
    print(f"Program has been running for {elapsed_time:.1f} seconds")
    time.sleep(0.1)
