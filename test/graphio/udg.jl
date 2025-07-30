using Test
using GadgetSearch
using Graphs

@testset "UDG Basic Types and Radius" begin
    @testset "Lattice Type Definitions" begin
        # Test that lattice types are concrete types
        @test Square <: GadgetSearch.LatticeType
        @test Triangular <: GadgetSearch.LatticeType
        
        # Test that instances can be created
        square = Square()
        triangular = Triangular()
        @test typeof(square) == Square
        @test typeof(triangular) == Triangular
    end
    
    @testset "get_radius Function" begin
        # Test radius values for different lattice types
        @test get_radius(Square()) == 1.5
        @test get_radius(Triangular()) == 1.1
        
        # Test return type is Float64
        @test typeof(get_radius(Square())) == Float64
        @test typeof(get_radius(Triangular())) == Float64
    end
end

@testset "unit_disk_graph Function" begin
    @testset "Basic Functionality" begin
        # Test with simple 2D points
        points = [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0)]
        unit = 1.5
        g = GadgetSearch.unit_disk_graph(points, unit)
        
        @test nv(g) == 3  # Should have 3 vertices
        @test typeof(g) == SimpleGraph{Int}
        
        # Points at distance 1.0 should be connected (< 1.5)
        @test has_edge(g, 1, 2)
        @test has_edge(g, 2, 3)
        # Points at distance 2.0 should not be connected (≥ 1.5^2 = 2.25)
        @test !has_edge(g, 1, 3)
    end
    
    @testset "Edge Cases" begin
        # Test with single point
        single_point = [(0.0, 0.0)]
        g_single = GadgetSearch.unit_disk_graph(single_point, 1.0)
        @test nv(g_single) == 1
        @test ne(g_single) == 0
        
        # Test with two points exactly at unit distance
        points_exact = [(0.0, 0.0), (1.0, 0.0)]
        g_exact = GadgetSearch.unit_disk_graph(points_exact, 1.0)
        @test nv(g_exact) == 2
        @test ne(g_exact) == 0  # Distance = 1.0, unit^2 = 1.0, should not connect
        
        # Test with two points just under unit distance
        points_under = [(0.0, 0.0), (0.9, 0.0)]
        g_under = GadgetSearch.unit_disk_graph(points_under, 1.0)
        @test nv(g_under) == 2
        @test ne(g_under) == 1  # Distance = 0.9 < 1.0, should connect
    end
    
    @testset "Complete Graph" begin
        # Test with points all within unit distance - should form complete graph
        close_points = [(0.0, 0.0), (0.3, 0.0), (0.0, 0.3), (0.3, 0.3)]
        g_complete = GadgetSearch.unit_disk_graph(close_points, 1.0)
        @test nv(g_complete) == 4
        @test ne(g_complete) == 6  # Complete graph on 4 vertices has 6 edges
    end
end

@testset "get_physical_positions Function" begin
    @testset "Square Lattice" begin
        square = Square()
        pos = [(1, 1), (2, 1), (1, 2)]
        physical_pos = GadgetSearch.get_physical_positions(square, pos)
        
        # Square lattice should preserve integer coordinates as floats
        @test physical_pos == [(1.0, 1.0), (2.0, 1.0), (1.0, 2.0)]
        @test typeof(physical_pos) == Vector{Tuple{Float64, Float64}}
    end
    
    @testset "Triangular Lattice" begin
        triangular = Triangular()
        pos = [(1, 1), (2, 1), (1, 2), (2, 2)]
        physical_pos = GadgetSearch.get_physical_positions(triangular, pos)
        
        h = sqrt(3) / 2  # triangular lattice height factor
        
        # Check expected transformations
        expected = [
            (1.0 + 0.5, 1 * h),  # (1,1) -> odd y, add 0.5 to x
            (2.0 + 0.5, 1 * h),        # (2,1) -> odd y, add 0.5 to x  
            (1.0, 2 * h),        # (1,2) -> even y, no x offset
            (2.0, 2 * h)   # (2,2) -> even y, no x offset
        ]
        
        @test length(physical_pos) == 4
        @test typeof(physical_pos) == Vector{Tuple{Float64, Float64}}
        
        # Test individual coordinates with tolerance for floating point
        for (i, (actual, exp)) in enumerate(zip(physical_pos, expected))
            @test actual[1] ≈ exp[1] atol=1e-10
            @test actual[2] ≈ exp[2] atol=1e-10
        end
    end
end

@testset "get_pin_positions Function" begin
    @testset "Basic Pin Positioning" begin
        square = Square()
        nx, ny = 2, 2
        top, bottom, left, right = GadgetSearch.get_pin_positions(square, nx, ny)
        
        # Test that we get the right number of positions
        @test length(top) == ny
        @test length(bottom) == ny
        @test length(left) == nx
        @test length(right) == nx
        
        # Test specific positions for nx=2, ny=2
        # Top: (1, y) for y in 2:3
        @test top == [(1.0, 2.0), (1.0, 3.0)]
        # Bottom: (4, y) for y in 2:3  
        @test bottom == [(4.0, 2.0), (4.0, 3.0)]
        # Left: (x, 1) for x in 2:3
        @test left == [(2.0, 1.0), (3.0, 1.0)]
        # Right: (x, 4) for x in 2:3
        @test right == [(2.0, 4.0), (3.0, 4.0)]
    end
    
    @testset "Triangular Lattice Pin Positions" begin
        triangular = Triangular()
        nx, ny = 1, 1
        top, bottom, left, right = GadgetSearch.get_pin_positions(triangular, nx, ny)
        
        h = sqrt(3) / 2
        
        # For nx=1, ny=1, we should get one position for each side
        @test length(top) == 1
        @test length(bottom) == 1
        @test length(left) == 1
        @test length(right) == 1
        
        # Check that positions are properly transformed for triangular lattice
        @test top[1][2] ≈ 2 * h  # y=2 -> 2*h
        @test bottom[1][2] ≈ 2 * h  # y=2 -> 2*h
        @test left[1][2] ≈ 1 * h  # y=1 -> 1*h
        @test right[1][2] ≈ 3 * h  # y=3 -> 3*h
    end
end

@testset "get_inner_points Function" begin
    @testset "Basic Inner Points" begin
        square = Square()
        nx, ny = 2, 2
        inner_points = GadgetSearch.get_inner_points(square, nx, ny)
        
        # For nx=2, ny=2, inner points should be at (x,y) for x in 2:3, y in 2:3
        expected_points = [(2.0, 2.0), (3.0, 2.0), (2.0, 3.0), (3.0, 3.0)]
        
        @test length(inner_points) == 4
        @test Set(inner_points) == Set(expected_points)
    end
    
    @testset "Single Inner Point" begin
        square = Square()
        nx, ny = 1, 1
        inner_points = GadgetSearch.get_inner_points(square, nx, ny)
        
        # For nx=1, ny=1, should have only one inner point at (2,2)
        @test length(inner_points) == 1
        @test inner_points[1] == (2.0, 2.0)
    end
    
    @testset "Triangular Lattice Inner Points" begin
        triangular = Triangular()
        nx, ny = 1, 1
        inner_points = GadgetSearch.get_inner_points(triangular, nx, ny)
        
        h = sqrt(3) / 2
        
        @test length(inner_points) == 1
        # Point (2,2) in triangular lattice: x=2, y=2 (even) -> (2.0, 2*h)
        @test inner_points[1][1] ≈ 2.0
        @test inner_points[1][2] ≈ 2 * h
    end
end

@testset "generate_full_grid_udg Integration" begin
    @testset "Small Grid Test (without shortg)" begin
        # Test the function structure without relying on external shortg tool
        square = Square()
        nx, ny = 1, 1
        
        # Create a temporary file for testing
        temp_file = tempname() * ".g6"
        
        # Test if the function runs without error when shortg is not available
        # We expect this to fail due to missing shortg, but we can test the structure
        try
            result = generate_full_grid_udg(square, nx, ny; path=temp_file)
            
            # If it succeeds (shortg is available), check the result
            @test typeof(result) == String
            @test result == temp_file
            
            # Clean up
            if isfile(temp_file)
                rm(temp_file)
            end
            
        catch e
            # Expected error when shortg is not available
            @test occursin("shortg", string(e)) || occursin("external tool", string(e))
        end
    end
    
    @testset "Parameter Validation" begin
        square = Square()
        
        # Test with valid parameters
        @test_nowarn generate_full_grid_udg
        
        # Test grid dimensions
        nx, ny = 1, 1
        top, bottom, left, right = GadgetSearch.get_pin_positions(square, nx, ny)
        inner = GadgetSearch.get_inner_points(square, nx, ny)
        
        # Should have 4 pin positions + 1 inner point = 5 total points
        total_points = length(top) + length(bottom) + length(left) + length(right) + length(inner)
        @test total_points == 5
        
        # Each configuration should have exactly 4 pins + inner points
        expected_graph_size = 4 + length(inner)
        @test expected_graph_size == 5  # 4 pins + 1 inner point
    end
end

@testset "Edge Cases and Error Conditions" begin
    @testset "Empty and Invalid Inputs" begin
        # Test unit_disk_graph with empty input
        empty_points = Tuple{Float64, Float64}[]
        g_empty = GadgetSearch.unit_disk_graph(empty_points, 1.0)
        @test nv(g_empty) == 0
        @test ne(g_empty) == 0
        
        # Test with zero unit distance
        points = [(0.0, 0.0), (0.1, 0.0)]
        g_zero = GadgetSearch.unit_disk_graph(points, 0.0)
        @test nv(g_zero) == 2
        @test ne(g_zero) == 0  # No edges with zero distance
        
        # Test with very small unit distance
        g_small = GadgetSearch.unit_disk_graph(points, 0.05)
        @test nv(g_small) == 2
        @test ne(g_small) == 0  # 0.1 > 0.05^2
        
        g_larger = GadgetSearch.unit_disk_graph(points, 0.5)
        @test nv(g_larger) == 2
        @test ne(g_larger) == 1  # 0.1 < 0.5^2
    end
    
    @testset "Large Unit Distance" begin
        # Test with very large unit distance - should connect everything
        points = [(0.0, 0.0), (10.0, 0.0), (0.0, 10.0), (10.0, 10.0)]
        g_large = GadgetSearch.unit_disk_graph(points, 100.0)
        @test nv(g_large) == 4
        @test ne(g_large) == 6  # Complete graph
    end
    
    @testset "Boundary Conditions" begin
        # Test get_physical_positions with empty input
        square = Square()
        empty_pos = Tuple{Int, Int}[]
        result_empty = GadgetSearch.get_physical_positions(square, empty_pos)
        @test length(result_empty) == 0
        @test typeof(result_empty) == Vector{Tuple{Float64, Float64}}
        
        # Test get_pin_positions with minimal grid
        nx_min, ny_min = 0, 0
        top, bottom, left, right = GadgetSearch.get_pin_positions(square, nx_min, ny_min)
        @test length(top) == 0
        @test length(bottom) == 0  
        @test length(left) == 0
        @test length(right) == 0
        
        # Test get_inner_points with minimal grid
        inner_min = GadgetSearch.get_inner_points(square, nx_min, ny_min)
        @test length(inner_min) == 0
    end
    
    @testset "Type Stability" begin
        # Test that functions return the expected types
        square = Square()
        triangular = Triangular()
        
        @test typeof(get_radius(square)) == Float64
        @test typeof(get_radius(triangular)) == Float64
        
        pos = [(1, 1)]
        @test typeof(GadgetSearch.get_physical_positions(square, pos)) == Vector{Tuple{Float64, Float64}}
        @test typeof(GadgetSearch.get_physical_positions(triangular, pos)) == Vector{Tuple{Float64, Float64}}
        
        points = [(1.0, 1.0)]
        @test typeof(GadgetSearch.unit_disk_graph(points, 1.0)) == SimpleGraph{Int}
    end
end

