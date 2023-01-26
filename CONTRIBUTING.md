# Contributing to `volesti`

:+1::tada: First off, thanks for taking the time to contribute! :tada::+1:

The following is a set of guidelines for contributing to volesti,
which are hosted in the [GeomScale Organization](https://github.com/GeomScale) on GitHub.
These are mostly guidelines, not rules.
Use your best judgment, and feel free to propose changes to this document in a pull request.

## Table of Contents

- [Prerequisites (how to start)](#prerequisites--how-to-start-)
- [Testing the development branch of volesti (get the tools ready)](#testing-the-development-branch-of-volesti--get-the-tools-ready-)
- [Fork volesti repository (this is your repo now!)](#fork--volesti-repository--this-is-your-repo-now--)
  - [Verify if your fork works (optional)](#verify-if-your-fork-works--optional-)
- [Working with volesti (get ready to contribute)](#working-with-volesti--get-ready-to-contribute-)
  - [GitFlow workflow](#gitflow-workflow)
  - [Create new branch for your work](#create-new-branch-for-your-work)
  - [Verify your new branch (optional)](#verify-your-new-branch--optional-)
- [Modify the branch (implement, implement, implement)](#modify-the-branch--implement--implement--implement-)
  - [Tests](#tests)
  - [Push](#push)
- [Pull request (the joy of sharing)](#pull-request--the-joy-of-sharing-)
- [Review (ok this is not an exam)](#review--ok-this-is-not-an-exam-)

## Prerequisites (how to start)

- git (see [Getting Started with Git](https://help.github.com/en/github/using-git/getting-started-with-git-and-github))
- a compiler to run tests - gcc, clang, etc.
- configured GitHub account

Other helpful links:

- http://git-scm.com/documentation
- https://help.github.com/articles/set-up-git
- https://opensource.com/article/18/1/step-step-guide-git

## Testing the development branch of volesti (get the tools ready)

Clone the repository,

    git clone git@github.com:GeomScale/volume_approximation.git volesti
    cd volesti
    git branch -vv

the last command should tell you that you are in `develop` branch.

To compile the `C++` code you have to specify the path to external library `liblpsolve55.so/dll/dylib` (see [here](doc/cpp_interface.md) more detail), by running, in folder test:

    mkdir -p test/build && cd test/build
    cmake -DLP_SOLVE=_PATH_TO_LIB_FILE_ ..
    # e.g. on linux: cmake -DLP_SOLVE=/usr/lib/lp_solve/liblpsolve55.so ..
    make

Run the tests,

    ctest -jK

where `K` is the number of CPU threads. By adding the option `--verbose` to `ctest` you get more information about the tests,
_e.g._ time per test, volume computed and the name of the polytope or convex body.

![test_cube](https://user-images.githubusercontent.com/3660366/72348403-0524df00-36e3-11ea-9b6d-288a2bddc22c.png)

If everything works for you, you may move forward.

## Fork volesti repository (this is your repo now!)

You can't work directly in the original volesti repository, therefore you should create your fork of this library.
This way you can modify the code and when the job is done send a pull request to merge your changes with the original
repository.

![fork](https://user-images.githubusercontent.com/3660366/72348562-57fe9680-36e3-11ea-9746-385ff61c752a.png)

1. login on `GitHub`
2. go to [volesti repository](https://github.com/GeomScale/volume_approximation)
3. click the 'Fork' button
4. choose your profile
5. wait
6. ready to contribute!

More info: [Forking Projects](https://guides.github.com/activities/forking/)

### Verify if your fork works (optional)

Go out of `volesti` directory

    cd ..

clone your repository and checkout develop branch

    git clone git@github.com:vissarion/volume_approximation.git volesti_fork
    cd volesti_fork
    git remote add upstream git@github.com:GeomScale/volesti.git
    git fetch upstream
    git checkout upstream/develop
    git branch -vv

see commits

    git log
    gitk

For now you should see exactly the same commits as in `volesti` repository.

## Working with volesti (get ready to contribute)

### GitFlow workflow

Volesit is using the [GitFlow](http://nvie.com/posts/a-successful-git-branching-model/) workflow.
It's because it is very well suited to collaboration and scaling the development team.
Each repository using this model should contain two main branches:

- master - release-ready version of the library
- develop - development version of the library

and could contain various supporting branches for new features and hotfixes.

As a contributor you'll most likely be adding new features or fixing bugs in the development version of the library.
This means that for each contribution you should create a new branch originating from the develop branch,
modify it and send a pull request in order to merge it, again with the develop branch.

### Create new branch for your work

Make sure you're in develop branch running

    git branch -vv

you should see

![branch -vv](https://user-images.githubusercontent.com/3660366/72348696-a1e77c80-36e3-11ea-93ec-70f5622c0675.png)

Now you should pick a name for your new branch that doesn't already exist.
The following checks for existing remote branches

    git branch -a

![List of branches](https://user-images.githubusercontent.com/3660366/72348763-c5aac280-36e3-11ea-8f2c-c66e2c107929.png)
Alternatively, you can check them on `GitHub`.

Assume you want to add some new functionality (i.e. a new feature) for example a new sampling algorithm. Then you have
to create a new branch e.g. `feature/the_fastest_sampling_algo_ever`

Create new local branch

    git branch feature/the_fastest_sampling_algo_ever
    git checkout feature/the_fastest_sampling_algo_ever

push it to your fork

    git push -u my_fork feature/the_fastest_sampling_algo_ever

Note that the `-u` switch also sets up the tracking of the remote branch. Your new branch now is created!

### Verify your new branch (optional)

Now with the command

    git branch -vv

you see

![branch-picked](https://user-images.githubusercontent.com/3660366/72348881-09053100-36e4-11ea-8187-c5a2fc7004b2.png)

Note that without the `-u` switch you wouldn't see the tracking information for your new branch.

Alternatively, your newly created remote branch is also available on GitHub

![new-feature-branch-github](https://user-images.githubusercontent.com/3660366/72349060-76b15d00-36e4-11ea-8065-e2367d5a2696.png)

## Modify the branch (implement, implement, implement)

Before contributiong to a library by adding a new feature, or a bugfix, or improving documentation,
it is always wise to interact with the community of developers, for example by opening an issue.

### Tests

Tests are placed in the `test` directory and use the [doctest](https://github.com/onqtam/doctest) library.

It is recommended to add new test whenever you contribute a new functionality/feature.
Also if your contribution is a bugfix then consider adding this case to the test-suite.

### Push

At the end, push your changes to the remote branch

    git push my_fork feature/the_fastest_sampling_algo_ever

or if your local branch is tracking the remote one, just

    git push

## Pull request (the joy of sharing)

After pushing your work you should be able to see it on `GitHub`.

Click "Compare and pull request" button or the "New pull request" button.

Add title and description

![RP](https://user-images.githubusercontent.com/3660366/72349397-21298000-36e5-11ea-9932-c8759c34ab2f.png)

and click the "Create pull request" button.

## Review (ok this is not an exam)

After creating a pull request your code will be reviewed. You can propose one or more reviewers
by clicking on the "Reviewers" button

![reviewer](https://user-images.githubusercontent.com/3660366/72349476-44ecc600-36e5-11ea-81cd-d0938d923529.png)

If there are no objections your changes will be merged.
Otherwise you'll see some comments under the pull request and/or under specific lines of your code.
Then you have to make the required changes, commit them and push to your branch.
Those changes will automatically be a part of the same pull request. This procedure will be repeated until the code
is ready for merging.

If you're curious how it looks like you may see one of the open or closed
[pull requests](https://github.com/GeomScale/volume_approximation/pulls).
