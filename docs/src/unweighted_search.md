# Dynamic unweighted search
`search_unweighted_gadgets` builds planar logical graphs, applies exact vertex
splits and even subdivisions, and embeds them on KSG or triangular lattices.
`is_gadget_replacement` remains final; four-pin results also pass G1--G4.
`result.trace` stores edge lists and rewrite actions, never graph6 identifiers.
