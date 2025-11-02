# Instructions for Creating Pull Request

## Step 1: Fork the Repository (if not already done)

1. Go to https://github.com/GeomScale/volesti
2. Click the "Fork" button in the top right
3. This creates a copy at https://github.com/YOUR_USERNAME/volesti

## Step 2: Add Your Fork as a Remote

Replace `YOUR_USERNAME` with your actual GitHub username:

```bash
git remote add fork https://github.com/YOUR_USERNAME/volesti.git
```

## Step 3: Push Your Branch to Your Fork

```bash
git push -u fork interface-to-GNU-Octave-language
```

## Step 4: Create Pull Request via GitHub Web Interface

1. Go to https://github.com/YOUR_USERNAME/volesti
2. You should see a banner saying "interface-to-GNU-Octave-language had recent pushes"
3. Click the "Compare & pull request" button
4. Or manually:
   - Go to https://github.com/GeomScale/volesti
   - Click "Pull requests" tab
   - Click "New pull request"
   - Select "compare across forks"
   - Base repository: `GeomScale/volesti`, base: `develop`
   - Head repository: `YOUR_USERNAME/volesti`, compare: `interface-to-GNU-Octave-language`
   - Click "Create pull request"

## Step 5: Fill in Pull Request Details

Use this title:
```
Add GNU Octave interface package for VolEsti
```

And this description (see `octave/PULL_REQUEST.md` for full details):

```markdown
## Summary

This PR adds a complete GNU Octave package interface for VolEsti, allowing Octave users to use the library with native Octave code.

## Features

- Volume computation with multiple algorithms
- Uniform sampling with various random walk methods
- Polytope generation (hypercubes)
- Comprehensive test suite
- Complete documentation

## Implementation

- MEX interfaces for C++ library (`volesti_volume_mex.cpp`, `volesti_sample_mex.cpp`)
- Octave wrapper functions (`volume()`, `sample_points()`, `GenCube()`)
- Test suite with 11+ test cases
- Build system with CMake
- Installation and usage documentation

See `octave/README.md` and `octave/PULL_REQUEST.md` for more details.
```

## Alternative: Using GitHub CLI (if installed)

If you have GitHub CLI (`gh`) installed:

```bash
# Authenticate (if not already done)
gh auth login

# Push to fork
git push -u fork interface-to-GNU-Octave-language

# Create PR
gh pr create --base develop --title "Add GNU Octave interface package for VolEsti" --body-file octave/PULL_REQUEST.md
```

