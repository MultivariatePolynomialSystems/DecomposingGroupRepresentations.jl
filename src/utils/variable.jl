struct Variable
    var::Expression

    function Variable(var::Expression)
        fs = free_symbols(var)
        @assert length(fs) == 1 && fs[1] == var
        return new(var)
    end
end

Base.show(io::IO, v::Variable) = show(io, v.var)