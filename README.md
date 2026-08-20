## Phobos-Earth Transit Predictor

A MATLAB project that predicts and visualizes when Phobos (Mars's inner moon) passes in front of the Sun as seen from Earth, built from scratch using real orbital mechanics.

### Overview

This project figures out the positions of Earth, Mars, and Phobos using real orbital data, finds when a transit will happen, and creates a sky-map showing what it would look like — without using any pre-built astronomy toolbox.

### Key Features

* Tracks the orbits of Mars and Phobos using real JPL orbital data, solving for their positions with a Newton-Raphson method.
* Uses the IAU Mars pole model to correctly line up Phobos's orbit with how it would appear from Earth.
* Searches for transit windows in two steps — a broad scan first, then a finer scan to nail down the exact timing without checking every possible moment at full detail.
* Produces a sky-map showing the transit as it would appear from Earth.

### Architecture

The project is split into 8 MATLAB files, each handling one part of the process — orbital data, coordinate conversions, the solver, the search, and the plotting — so each part can be tested and changed on its own.

### Motivation

I built this on my own to really understand orbital mechanics — coding the actual math myself (orbit tracking, coordinate conversions, solving for positions) instead of relying on a toolbox to do it for me.
