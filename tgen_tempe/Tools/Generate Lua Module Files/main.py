import json
import os

# Open module json
with open("module_info.json", "r") as f:
    modules = json.load(f)

# Loop through every module
for mod in modules:

    # Create file name and dir
    app_dir = mod["AppName"]
    file = f"{dir}/{mod['AppVersion']}-Container.lua"

    # Create dir if it does not already exist
    if not os.path.exists(app_dir):
        os.mkdir(app_dir)

    # Set variables that will not be changed
    load = 'load("singluarity")'
    execute = 'execute {cmd="shopt -s expand_aliases",modeA={"load"}}'

    # Setting up documentation variables
    docs = " \n".join(mod["help_docs"])
    docs = f"help([[\nFor documentation visit:\n{docs}\n]])"

    # Generating whatis commands
    whatiscmds = [f'whatis("{cmd}")' for cmd in mod["whatis"]]
    whatiscmds = " \n".join(whatiscmds)

    # Generating Description help command
    desc = " \n".join(mod["help_desc"])
    desc = f"help([==[\n{desc}\n]==])"

    # Generating set_alias commands
    aliascmds = [f'set_alias("{k}","{v}")' for k, v in mod["cmd_aliases"].items()]
    aliases = " \n".join(aliascmds)

    # Generate data for file
    data = f"{docs}\n\n{whatiscmds}\n\n{desc}\n\n{load}\n{execute}\n{aliases}"

    # Create file and write data
    with open(file, "w+") as f:
        f.write(data)