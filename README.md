# PedestFlow

## About the model

This code tries to simulate a pedestrian flow in two crossing streets, as explained [here](https://gist.github.com/meiemari/b20406e05eb2aed6367361f85d552802​​​).

## Equations of motion

The rules of the pedestrian behaviour are encoded into the following equations of motion:

$$
\frac{dr}{dt} = \boldsymbol{v}_i(t)
$$

## Running locally

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> PedestFlow

It is authored by Ehsan Irani.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
