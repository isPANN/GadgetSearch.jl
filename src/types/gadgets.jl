"""Abstract type for all gadget types in the system"""
abstract type AbstractGadget end

"""Represents a gadget with its properties"""
struct Gadget{T<:Real} <: AbstractGadget
    rule_id::Int
    ground_states::Vector{String}
    graph_id::Int
    graph::SimpleGraph{Int}
    pins::Vector{Int}
    weights::Vector{T}
end 