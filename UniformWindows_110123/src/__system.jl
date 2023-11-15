""" Serializable representation of a physical system. """
module Systems
    import ..Float
        #= NOTE: This can easily be replaced by `const Float = Float64`
                to make this module totally independent of a job. =#

    import LinearAlgebra
    import Serialization
    import NPZ

    const MODEL_PATH = "$(ENV["HOME"])/CtrlJobs/matrix/model"
    const SECTORPATH = "$(ENV["HOME"])/CtrlJobs/matrix/sector"
    const EIGEN_PATH = "$(ENV["HOME"])/CtrlJobs/matrix/eigen"


    """
        Model(args...)

    Encapsulate a problem Hamiltonian, its diagonalization, and quantum number matrices.

    """
    struct Model
        H::Matrix{Complex{Float}}
        n::Int
        Λ::Vector{Float}
        U::Matrix{Complex{Float}}
        N::Matrix{Complex{Float}}
        S2::Matrix{Complex{Float}}
        Sz::Matrix{Complex{Float}}
    end

    """
        System(model, ψ0)

    Encapsulate a Model and a reference statevector, with useful scalar attributes.

    """
    struct System
        model::Model
        ψ0::Vector{Complex{Float}}

        # FOR CONVENIENT ACCESS
        n::Int
        REF::Float
        FCI::Float

        System(model::Model, ψ0::Vector{Complex{Float}}) = new(
            model,
            ψ0,
            # FOR CONVENIENT ACCESS
            model.n,
            real(LinearAlgebra.dot(ψ0, model.H*ψ0)),
            model.Λ[1],
        )
    end










    function Model(code::String; load_ΛU=true, save_ΛU=true)
        modelfile = "$MODEL_PATH/$code.npy"
        eigenfile = "$EIGEN_PATH/$code.ΛU"
        eigenfile_exists = isfile(eigenfile)

        suffix = get_suffix(code)
        sectorfile(prefix) = "$SECTORPATH/$(prefix)_$(suffix).npy"

        # LOAD MATRIX
        H = NPZ.npzread(modelfile)

        # INFER NUMBER OF QUBITS
        n = round(Int, log2(size(H,1)))

        # LOAD OR CONSTRUCT DIAGONALIZATION OF H
        ΛU = load_ΛU && eigenfile_exists ?
                Serialization.deserialize(eigenfile) :
                LinearAlgebra.eigen(LinearAlgebra.Hermitian(H))

        # SAVE DIAGONALIZATION, IF IT DOESN'T ALREADY EXIST
        if save_ΛU && !(load_ΛU && eigenfile_exists)
            Serialization.serialize(eigenfile, ΛU)
        end

        # UNPACK DIAGONALIZATION INTO EIGENVALUES AND EIGENVECTORS
        Λ, U = ΛU

        # FETCH SECTOR AND QUANTUM NUMBER OPERATORS
        N  = NPZ.npzread(sectorfile("N_"))
        S2 = NPZ.npzread(sectorfile("Sz"))
        Sz = NPZ.npzread(sectorfile("S2"))

        return Model(H, n, Λ, U, N, S2, Sz)
    end

    function System(code::String; load_ΛU=true, save_ΛU=true)
        prefix = get_prefix(code)
        suffix = get_suffix(code)
        ref_file = "$SECTORPATH/$(prefix)_$(suffix).npy"
        return System(
            Model(code; load_ΛU=load_ΛU, save_ΛU=save_ΛU),
            NPZ.npzread(ref_file),
        )
    end

    function get_prefix(code::String)
        for (pattern, specifier) in prefix_specifiers
            regex_match = match(pattern, code)
            isnothing(regex_match) && continue
            return specifier(regex_match)
        end
        error("Code pattern not recognized.")
    end

    get_suffix(L, N, M, mapping) = "$(L)_$(N)_$(M)_$(mapping)"

    function get_suffix(code::String)
        for (pattern, specifier) in suffix_specifiers
            regex_match = match(pattern, code)
            isnothing(regex_match) && continue
            return specifier(regex_match)
        end
        error("Code pattern not recognized.")
    end


    ######################################################################################
    #= The remaining mess chooses the sector file names by decoding the model codes. =#

    """ I have typically been using a standard abbreviation for mappings,
            but my model codes abbreviate them even more. """
    function code_to_mapping(map, tap)
        mapping = (
            map == "J" ? "JW" :
            map == "B" ? "BK" :
            map == "P" ? "P" :
            error("Unrecognized mapping code $map")
        )
        isnothing(tap) || (mapping *= "-$tap")
        return mapping
    end

    """ Map regex patterns of a model code onto a function,
            which inputs a RegexMatch object
            and outputs a suffix for the sector. """
    const suffix_specifiers = Dict{Regex, Function}(
        # HYDROGEN
        r"^(c?H)(\d+)(.)(P|J|B)(m|n)?$" => (regex_match) -> (
            (cfg, num, geo, map, tap) = regex_match;

            N = parse(Int, num);
            L = N << 1;
            M = (N >> 1) + (N & 1);     # pyscf defaults to odd electrons with α-spin
            mapping = code_to_mapping(map, tap);

            get_suffix(L, N, M, mapping)
        ),
        # HUBBARD
        r"^(c?L)(\d+)(.)(.)(P|J|B)(m|n)?$" => (regex_match) -> (
            (cfg, num, geo, bas, map, tap) = regex_match;

            N = parse(Int, num);
            L = N << 1;
            M = N >> 1;                 # I default to odd electrons with β-spin
            mapping = code_to_mapping(map, tap);

            get_suffix(L, N, M, mapping)
        ),
        # OTHERS
        r"lih30" => (regex_match) -> get_suffix(6, 4, 2, "P-m"),
            #= NOTE: These are not really the correct mapping!
                The "lih30.npy" matrix was generated using qiskit's tapering,
                    not openfermion's, so these sector matrices are totally wrong.
                    Ignore them when studying this molecule! =#
        r"H215" => (regex_match) -> get_suffix(4, 2, 1, "P-m"),
    )

    """ Map regex patterns of a model code onto a function,
            which inputs a RegexMatch object
            and outputs a prefix for the reference state. """
    const prefix_specifiers = Dict{Regex, Function}(
        # HYDROGEN
        r"^(c?H)(\d+)(.)(P|J|B)(m|n)?$" => (regex_match) -> "REF",
        # HUBBARD
        r"^(c?L)(\d+)(.)(.)(P|J|B)(m|n)?$" => (regex_match) -> (
            (cfg, num, geo, bas, map, tap) = regex_match;
            bas == "A" ? "ANTI" : "REF";
        ),
        # OTHERS
        r"lih30" => (regex_match) -> "QKA",
        r"H215" => (regex_match) -> "REF",
    )

end