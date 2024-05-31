#= Set the standard styles for each dataset we'll be plotting. =#

module lih30_surveys
    export nodal, optimal, bisectal, uniform
    export harmonics, harmonicz, harmonicxyz

    using AdaptFigures_052324

    MODALJOBS = "../AdaptiveModes_041724/jobs"
    SQUAREJOBS = "../AdaptiveWindows_040624/jobs"

    #= SURVEY TUPLES =#

    sub_naive = (
        job = SQUARE,
        surveydir = "$SQUAREJOBS/lih30_uniform.each",
        eMET = 21.0,
        label = "Div/Naive",
        style = (
            color=2,
            shape=:square,
        ),
    )

    sub_adapt = (
        job = SQUARE,
        surveydir = "$SQUAREJOBS/lih30_bisection.one",
        eMET = 24.0,
        label = "Div/Adaptive",
        style = (
            color=1,
            shape=:x,
        ),
    )

    sub_optim = (   # This is a bad abbreviation... :)
        job = SQUARE,
        surveydir = "$SQUAREJOBS/lih30_nodes.one",
        eMET = 22.0,
        label = "Div/Optimal",
        style = (
            color=3,
            shape=:utriangle,
        ),
    )

    sup_naive = (
        job = MODAL,
        surveydir = "$MODALJOBS/lih30_harmonics.complex.iterative",
        eMET = 24.0,
        label = "Sup/Naive",
        style = (
            color=2,
            shape=:circle,
        ),
    )

    sup_adapt = (
        job = MODAL,
        surveydir = "$MODALJOBS/lih30_harmonics",
        eMET = 22.0,
        label = "Sup/Adaptive",
        style = (
            color=1,
            shape=:+,
        ),
    )

    sup_optim = (
        job = MODAL,
        surveydir = "$MODALJOBS/lih30_gradients.modalharmonics",
        eMET = 22.0,
        label = "Sup/Optimal",
        style = (
            color=3,
            shape=:dtriangle,
        ),
    )






    nodal = (
        job = SQUARE,
        surveydir = "$SQUAREJOBS/lih30_nodes.one",
        eMET = 22.0,
        label = "Optimal",
        style = (
            color=1,
            shape=:utriangle,
        ),
    )

    optimal = (
        job = SQUARE,
        surveydir = "$SQUAREJOBS/lih30_optimal",
        eMET = 22.0,
        label = "Optimal (Coupled)",
        style = (
            color=2,
            shape=:utriangle,
        ),
    )

    bisectal = (
        job = SQUARE,
        surveydir = "$SQUAREJOBS/lih30_bisection",
        eMET = 24.0,
        label = "Bisectal",
        style = (
            color=3,
            shape=:+,
        ),
    )

    uniform = (
        job = MODAL,
        surveydir = "$MODALJOBS/lih30_square.complex.iterative.frozen",
        eMET = 24.0,
        label = "Uniform (Iterative)",
        style = (
            color=4,
            shape=:square,
        ),
    )

    harmonics = (
        job = MODAL,
        surveydir = "$MODALJOBS/lih30_harmonics",
        eMET = 22.0,
        label = "Harmonic",
        style = (
            color=5,
            shape=:circle,
        ),
    )

    harmonicz = (   # z for conventional z = x + iy
        job = MODAL,
        surveydir = "$MODALJOBS/lih30_harmonics.complex",
        eMET = 23.0,
        label = "Harmonic (Coupled)",
        style = (
            color=6,
            shape=:circle,
        ),
    )

    harmonicxyz = (
        job = MODAL,
        surveydir = "$MODALJOBS/lih30_harmonics.complex.iterative",
        eMET = 24.0,
        label = "Harmonic (Iterative)",
        style = (
            color=7,
            shape=:diamond,
        ),
    )



    nodal_each = (
        job = SQUARE,
        surveydir = "$SQUAREJOBS/lih30_nodes.oneeach",
        eMET = 23.0,
        label = "Optimal (Each Pulse)",
        style = (
            color=8,
            shape=:ltriangle,
            markersize=7,
        ),
    )

    nodal_all = (
        job = SQUARE,
        surveydir = "$SQUAREJOBS/lih30_nodes.oneall",
        eMET = 21.0,
        label = "Optimal (Each Window)",
        style = (
            color=9,
            shape=:rtriangle,
            markersize=7,
        ),
    )




end

