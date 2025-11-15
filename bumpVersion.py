from sys import argv

files = ["CMakeLists.txt", "src/infrastructure/ui/ui.cpp", "README.md"]
    
def main():
    current_version = get_current_version()
    major = int(current_version[0])
    minor = int(current_version[1])
    patch = int(current_version[2])

    new_major, new_minor, new_patch = update_version(major, minor, patch)
    update_files(new_major, new_minor, new_patch, major, minor, patch)


def get_current_version():
    with open(files[0], 'r') as f:
        lines = f.readlines()
        project_start = -1
        project_end = -1
        for line in lines:
            if "project" in line:
                project_start = lines.index(line)
            if ")" in line and project_start != -1:
                project_end = lines.index(line)
            if project_start != -1 and project_end != -1:
                break

        for i in range(project_start, project_end + 1):
            if "VERSION" in lines[i]:
                args = lines[i].strip().split(" ")
                version = args[1].split(".")
                return version


def update_version(major, minor, patch):
    if len(argv) == 1:
        print('Usage: python3 bumpVersion.py [+-][major|minor|patch]')
        print(f'Current version: {major}.{minor}.{patch}')
    
    if len(argv) == 2:
        print(f'Current version: {major}.{minor}.{patch}')
        if argv[1] == '+major':
            major += 1
        elif argv[1] == '+minor':
            minor += 1
        elif argv[1] == '+patch':
            patch += 1
        elif argv[1] == '-major':
            assert(major > 0)
            major -= 1
        elif argv[1] == '-minor':
            assert(minor > 0)
            minor -= 1
        elif argv[1] == '-patch':
            assert(patch > 0)
            patch -= 1
        else:
            print('Usage: python3 bumpVersion.py [+-][major|minor|patch]')
            print(f'Current version: {major}.{minor}.{patch}')
            return
        
        print(f'New version: {major}.{minor}.{patch}')

    return major, minor, patch


def update_files(new_major, new_minor, new_patch, major, minor, patch):
    for file in files:
        with open(file, 'r') as f:
            lines = f.readlines()
            for i in range(0, len(lines)):
                if f'{major}.{minor}.{patch}' in lines[i]:
                    lines[i] = lines[i].replace(f'{major}.{minor}.{patch}', f'{new_major}.{new_minor}.{new_patch}')
        
        with open(file, 'w') as f:
            for line in lines:
                f.write(line)


if __name__ == "__main__":
    main()
