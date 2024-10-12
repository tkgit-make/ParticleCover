Branch containing the folder DarkMatterDirectDetection -- this points to the git submodule containing the code.

Since we are using git submodules, to clone use the `--recursive` option:
```
git clone --recursive <this repo>
```

If by accident you forget to use `--recursive` in cloning (will create an empty dir), or simply enter the directory and call the `submodule update` command:
```
cd <name of directory>
git submodule update --init --recursive
```

To update, run `git pull` inside the submodule:
```
cd <name of directory>
git pull
```

At time of writing `<name of directory>` is DarkMatterDirectDetection.
