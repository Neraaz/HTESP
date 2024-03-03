"""Written by Niraj K. Nepal, Ph.D."""
import sys
import os

def rewrite_script(script1: str, script2: str, func_name: str) -> None:
    """
    Rewrite a Bash script with a modified function and check its syntax.

    Args:
    - script1 (str): Path to the original Bash script.
    - script2 (str): Path to the new Bash script to be created.
    - func_name (str): Name of the function to be modified in the script.

    Returns:
    - None
    """
    # Read the content of the original script
    with open(script1, "r") as read_script:
        lines = read_script.readlines()

    # Find the indices of the function block
    func_index = []
    for i, line in enumerate(lines):
        if ('for ((' in line or 'done' in line) and "all" not in line:
            func_index.append(i)

    # Extract the header and function body
    header = lines[0:func_index[0]]
    lhead = len(header)
    file_dft = header[lhead - 2:lhead]
    header = header[:lhead - 2]
    func_string = lines[func_index[0] + 1:func_index[-1]]

    # Write the modified script
    with open(script2, "w") as write_new:
        for head in header:
            write_new.write(head)
        write_new.write(func_name + "() {\n")
        write_new.write("    # Assign parameters to local variables\n")
        write_new.write("    ii=$1\n")
        write_new.write("    dft=$2\n")
        write_new.write("    file=$3\n")
        for func in func_string:
            write_new.write("    " + func)
        write_new.write("}\n")
        write_new.write("# Main section of the script\n")
        write_new.write("""echo "---------------------------------------------------------------------------------------------------------------"\n""")
        write_new.write("""echo "Submitting crystal structure relaxation"\n""")
        write_new.write("""echo "---------------------------------------------------------------------------------------------------------------"\n""")
        for dft in file_dft:
            write_new.write(dft)
        write_new.write("# Check if parallel command is available\n")
        write_new.write("if command -v parallel &>/dev/null; then\n")
        write_new.write("  NUM_CORES=$(grep -c '^processor' /proc/cpuinfo)\n")
        write_new.write("  NUM_PROCESSES=$((NUM_CORES + 1))\n")
        write_new.write("\n")
        write_new.write("  # Export functions for parallel execution\n")
        write_new.write("  export -f vasprun\n")
        write_new.write("  export -f qerun\n")
        write_new.write(f"  export -f {func_name}\n")
        write_new.write("\n")
        write_new.write("  # Execute relaxation tasks in parallel\n")
        write_new.write("""  parallel -j "$NUM_PROCESSES" {} ::: $(seq "$1" "$2") ::: "$dft" ::: "$file"\n""".format(func_name))
        write_new.write("else\n")
        write_new.write("  # Parallel command not found, run relaxation tasks in serial mode\n")
        write_new.write("""  echo "Parallel command not found. Running in serial mode."\n""")
        write_new.write("""  for ii in $(seq "$1" "$2"); do\n""")
        write_new.write("""    {} "$ii" "$dft" "$file"\n""".format(func_name))
        write_new.write("  done\n")
        write_new.write("fi\n")
        write_new.write("""echo "all done"\n""")

    # Check the syntax of the generated script
    print(f"Checking {script2}\n")
    os.system(f"bash -n {script2}")
    print("Looking good\n")
    print("Making executable...\n")
    os.system(f"chmod +x {script2}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py script1 script2 func_name")
        sys.exit(1)

    script1 = sys.argv[1]
    script2 = sys.argv[2]
    func_name = sys.argv[3]

    rewrite_script(script1, script2, func_name)
