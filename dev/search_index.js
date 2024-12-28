var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = SlidingDFTs","category":"page"},{"location":"api/#Usage","page":"API","title":"Usage","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"sdft","category":"page"},{"location":"api/#SlidingDFTs.sdft","page":"API","title":"SlidingDFTs.sdft","text":"sdft(method, x[, safe=true])\n\nReturn an iterator to produce a sliding DFT of x using the given method. If safe == true (the default behavior), this iterator produces a new vector at each iteration.\n\nSet the optional argument safe=false to improve performance by reducing allocations, at the expense of unexpected behavior if the resulting vector is mutated between iterations.\n\n\n\n\n\n","category":"function"},{"location":"api/#Sliding-DFT-methods","page":"API","title":"Sliding DFT methods","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [SlidingDFTs]\nPages = [\"SDFT.jl\"]","category":"page"},{"location":"api/#SlidingDFTs.SDFT","page":"API","title":"SlidingDFTs.SDFT","text":"SDFT(n)\n\nBasic method to compute a Sliding Discrete Fourier Transform of window length n [1], through the recursive formula:\n\nX_i+1k = exp(j2pikn) cdot (X_ik + xi+n - xi)\n\nThe transfer function for the k-th bin of this method is:\n\n$ H(z) = \\frac{1 - z^{-n}}{1 - \\exp(j2{\\pi}k/n) z^{-1}} $\n\nReferences\n\n[1] Jacobsen, E. and Lyons, R. \"The sliding DFT,\" IEEE Signal Processing Magazine, 20(2), 74-80 (2003), DOI: 10.1109/MSP.2003.1184347\n\n\n\n\n\n","category":"type"},{"location":"api/#Implementation","page":"API","title":"Implementation","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [SlidingDFTs]\nPages = [\"interface.jl\"]\nFilter = !isequal(SlidingDFTs.sdft)","category":"page"},{"location":"api/#SlidingDFTs.dataoffsets-Tuple{Any}","page":"API","title":"SlidingDFTs.dataoffsets","text":"dataoffsets(method)\n\nReturn an integer or a vector of integers with the offsets of data samples that are needed by the given method to compute a sliding DFT.\n\nIf the code of SlidingDFTs.updatepdf! that dispatches on the type of method uses the function SlidingDFTs.previousdata, this function must return the integers that are used as the third argument (offset) of that function.\n\nIf that function is not needed (no past samples are used), this one may return nothing to reduce memory allocations.\n\n\n\n\n\n","category":"method"},{"location":"api/#SlidingDFTs.dftback-Tuple{Any}","page":"API","title":"SlidingDFTs.dftback","text":"dftback(method)\n\nReturn an integer or a vector of positive integers with the indices of the previous iterations that are needed by the given method to compute a sliding DFT.\n\nIf the code of SlidingDFTs.updatepdf! for the type of method uses the function SlidingDFTs.previousdft, this function must return the integers that are used as the third argument (back) of that function.\n\nIf that function is not needed, this one may return nothing to reduce memory allocations.\n\n\n\n\n\n","category":"method"},{"location":"api/#SlidingDFTs.iterationcount-Tuple{Any}","page":"API","title":"SlidingDFTs.iterationcount","text":"iterationcount(state)\n\nReturn the number of iterations done for the sliding DFT represented by state.\n\nIf the DFT computed in the most recent iteration corresponds to the fragment of the data series between its positions i and i+n, then this function returns  the number i.\n\n\n\n\n\n","category":"method"},{"location":"api/#SlidingDFTs.nextdata-Tuple{SlidingDFTs.StateData}","page":"API","title":"SlidingDFTs.nextdata","text":"nextdata(state)\n\nReturn the next value after the fragment of the data series that was used in the most recent iteration of the sliding DFT represented by state.\n\nIf the DFT computed in the most recent iteration corresponds to the fragment of the data series between its positions i and i+n, then this function returns the i+n+1-th value.\n\nThere is no defined behavior if such value does not exist (i.e. if the end of a finite data series was reached).\n\n\n\n\n\n","category":"method"},{"location":"api/#SlidingDFTs.previousdata","page":"API","title":"SlidingDFTs.previousdata","text":"previousdata(state[, offset=0])\n\nReturn the first value of the fragment of the data series that was used in the most recent iteration of the sliding DFT represented by state, or at offset positions after the beginning of that fragment.\n\nIf the DFT computed in the most recent iteration corresponds to the fragment of the data series between its positions i and i+n, then this function returns the i+offset-th value.\n\n\n\n\n\n","category":"function"},{"location":"api/#SlidingDFTs.previousdft","page":"API","title":"SlidingDFTs.previousdft","text":"previousdft(state[, back=0])\n\nReturn the DFT computed in the most recent iteration of the sliding DFT represented by state, or in a number of previous iterations equal to back.\n\nIf the DFT computed in the most recent iteration corresponds to the fragment of the data series between its positions i and i+n, then this function returns its DFT for the fragment between i-back and i+n-back.\n\n\n\n\n\n","category":"function"},{"location":"api/#SlidingDFTs.updatedft!","page":"API","title":"SlidingDFTs.updatedft!","text":"updatedft!(dft, x, method, state)\n\nUpdate the values of a sliding Discrete Fourier Transform (DFT) of a data series, according to the algorithm of the provided method, for a given state of the sliding DFT.\n\ndft is a mutable collection of complex values with length equal to windowlength(method), containing the value returned by the last iteration of the sliding DFT.\n\nx is the data series for which the sliding DFT is computed, at least as long as windowlength(method).\n\nmethod is the object that defines the method to compute the sliding DFT.\n\nstate is an object generated automatically at each iteration of an object created by sdft(method, x), containing the information that is needed to compute the new values of the sliding DFT, together with x. That information can be extracted from state with the following functions:\n\nSlidingDFTs.previousdft to get the DFTs of previous iterations.\nSlidingDFTs.previousdata to get a previous value of the data series.\nSlidingDFTs.nextdata to get the next value of the data series.\n\n\n\n\n\n","category":"function"},{"location":"api/#SlidingDFTs.windowlength","page":"API","title":"SlidingDFTs.windowlength","text":"SlidingDFTs.windowlength(method)\n\nReturn the length of the window used by method.\n\n\n\n\n\n","category":"function"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"CurrentModule = SlidingDFTs","category":"page"},{"location":"implement/#Impementing-SlidingDFTs","page":"Implementing SlidingDFTs","title":"Impementing SlidingDFTs","text":"","category":"section"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"SlidingDFTs provides an interface to define new methods to compute Sliding Discrete Fourier Transforms (SDFT), which can be implemented in pull requests to this package, but also in independent packages.","category":"page"},{"location":"implement/#Theoretical-basis","page":"Implementing SlidingDFTs","title":"Theoretical basis","text":"","category":"section"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"An SDFT is generally implemented as a recursive equation, such that if X_ik is the DFT of xj for j = i ldots i+n at frequency k, the next iteration is:","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"X_i+1k = f(k X_1k ldots X_ik xi ldots xi+n xi+n+1)","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"Such an equation depends on:","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"The frequency k\nThe values of the DFT in one or more previous iterations, X_pk for p = 1 ldots i.\nThe values of the data series used in the most recent iteration, xj for j = i ldots i+n.\nThe next value of the data series after the fragment used in the most recent iteration, xi+n+1.","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"For instance, the basic definition of the SDFT uses the following formula: ","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"X_i+1k = Wk cdot (X_ik + xi+n - xi)","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"which depends only on the most recent DFT (X_ik), the first data point of the fragment used in that DFT (xi), and the next data point after that fragment (xi+n+1), plus a \"twiddle\" factor Wk that only depends on the frequency k, equal to exp(j2pikn). Other variations of the SDFT may use formulas that depend on previous iterations or other values of the data series in that fragment.","category":"page"},{"location":"implement/#Implementation-in-Julia-object-types","page":"Implementing SlidingDFTs","title":"Implementation in Julia object types","text":"","category":"section"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"A method to compute an SDFT is defined by three kinds of object types:","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"One for the method, which contains the fixed parameters that are needed by the algorithm to compute the SDFT.\nAnother for the iterator created by the function sdft, which binds a method with the target data series.\nAnd yet another for the state of the iterator, which holds the information needed by the algorithm that depends on the data series and changes at each iteration.","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"The internals of SlidingDFTs take care of the design of the iterator and state types, and of the creation of their instances. The only thing that has to be defined to create a new kind of SDFT is the struct of the method with the fixed parameters, and a few function methods dispatching on that type. One of those functions is the one that implements the recursive formula, using also the information stored in the state.","category":"page"},{"location":"implement/#Extract-information-from-the-state-of-SDFT-iterators","page":"Implementing SlidingDFTs","title":"Extract information from the state of SDFT iterators","text":"","category":"section"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"Before explaining how to define a new type for SDFTs, it is convenient to know the functions that can be used to extract the information that is stored in the state of SDFT iterators. There is a function for each of the three kinds of variables used in the general recursive equation presented above.","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"previousdft\npreviousdata\nnextdata","category":"page"},{"location":"implement/#SlidingDFTs.previousdft-implement","page":"Implementing SlidingDFTs","title":"SlidingDFTs.previousdft","text":"previousdft(state[, back=0])\n\nReturn the DFT computed in the most recent iteration of the sliding DFT represented by state, or in a number of previous iterations equal to back.\n\nIf the DFT computed in the most recent iteration corresponds to the fragment of the data series between its positions i and i+n, then this function returns its DFT for the fragment between i-back and i+n-back.\n\n\n\n\n\n","category":"function"},{"location":"implement/#SlidingDFTs.previousdata-implement","page":"Implementing SlidingDFTs","title":"SlidingDFTs.previousdata","text":"previousdata(state[, offset=0])\n\nReturn the first value of the fragment of the data series that was used in the most recent iteration of the sliding DFT represented by state, or at offset positions after the beginning of that fragment.\n\nIf the DFT computed in the most recent iteration corresponds to the fragment of the data series between its positions i and i+n, then this function returns the i+offset-th value.\n\n\n\n\n\n","category":"function"},{"location":"implement/#SlidingDFTs.nextdata-implement","page":"Implementing SlidingDFTs","title":"SlidingDFTs.nextdata","text":"nextdata(state)\n\nReturn the next value after the fragment of the data series that was used in the most recent iteration of the sliding DFT represented by state.\n\nIf the DFT computed in the most recent iteration corresponds to the fragment of the data series between its positions i and i+n, then this function returns the i+n+1-th value.\n\nThere is no defined behavior if such value does not exist (i.e. if the end of a finite data series was reached).\n\n\n\n\n\n","category":"function"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"For instance, the values used in the formula of the basic SDFT may be obtained from a state object as:","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"SlidingDFTs.previousdft(state, 0) for X_i.\nSlidingDFTs.previousdata(state, 0) for xi.\nSlidingDFTs.nextdata(state) for xi+n+1.","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"Notice that the second arguments of previousdft and previousdata might have been ommited in this case, since they are zero by default.","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"For methods that need to know how many steps of the SDFT have been done, this can also be extracted:","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"iterationcount","category":"page"},{"location":"implement/#SlidingDFTs.iterationcount-implement","page":"Implementing SlidingDFTs","title":"SlidingDFTs.iterationcount","text":"iterationcount(state)\n\nReturn the number of iterations done for the sliding DFT represented by state.\n\nIf the DFT computed in the most recent iteration corresponds to the fragment of the data series between its positions i and i+n, then this function returns  the number i.\n\n\n\n\n\n","category":"function"},{"location":"implement/#Definition-of-new-SDFT-types","page":"Implementing SlidingDFTs","title":"Definition of new SDFT types","text":"","category":"section"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"The design of the struct representing a new SDFT type is free, but it is required to implement the following methods dispatching on that type:","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"windowlength\nupdatedft!","category":"page"},{"location":"implement/#SlidingDFTs.windowlength-implement","page":"Implementing SlidingDFTs","title":"SlidingDFTs.windowlength","text":"SlidingDFTs.windowlength(method)\n\nReturn the length of the window used by method.\n\n\n\n\n\n","category":"function"},{"location":"implement/#SlidingDFTs.updatedft!-implement","page":"Implementing SlidingDFTs","title":"SlidingDFTs.updatedft!","text":"updatedft!(dft, x, method, state)\n\nUpdate the values of a sliding Discrete Fourier Transform (DFT) of a data series, according to the algorithm of the provided method, for a given state of the sliding DFT.\n\ndft is a mutable collection of complex values with length equal to windowlength(method), containing the value returned by the last iteration of the sliding DFT.\n\nx is the data series for which the sliding DFT is computed, at least as long as windowlength(method).\n\nmethod is the object that defines the method to compute the sliding DFT.\n\nstate is an object generated automatically at each iteration of an object created by sdft(method, x), containing the information that is needed to compute the new values of the sliding DFT, together with x. That information can be extracted from state with the following functions:\n\nSlidingDFTs.previousdft to get the DFTs of previous iterations.\nSlidingDFTs.previousdata to get a previous value of the data series.\nSlidingDFTs.nextdata to get the next value of the data series.\n\n\n\n\n\n","category":"function"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"The formula of the basic SDFT formula could be implemented for a type MyBasicSDFT as follows:","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"import SlidingDFTs: updatedft!, windowlength, nextdata, previousdata\n\nfunction udpatedft!(dft, x, method::MyBasicSDFT, state)\n    n = windowlength(method)\n    for k in eachindex(dft)\n        X_i = dft[k]\n        x_iplusn = nextdata(state)\n        x_i = previousdata(state)\n        Wk = exp(2π*im*k/n)\n        dft[k] = Wk * (X_i + x_iplusn - x_i)\n    end\nend","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"(The type SDFT implemented in SlidingDFTs actually has as a similar, but not identic definition.)","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"Notice that it was not necessary to use previousdft to retrieve the most recent iteration of the DFT, since that is assummed to be already contained in the first argument dft.","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"Depending on the functions that are used in the particular implementation of updatepdf! for a given type, the following methods should be defined too:","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"dftback\ndataoffsets","category":"page"},{"location":"implement/#SlidingDFTs.dftback-implement","page":"Implementing SlidingDFTs","title":"SlidingDFTs.dftback","text":"dftback(method)\n\nReturn an integer or a vector of positive integers with the indices of the previous iterations that are needed by the given method to compute a sliding DFT.\n\nIf the code of SlidingDFTs.updatepdf! for the type of method uses the function SlidingDFTs.previousdft, this function must return the integers that are used as the third argument (back) of that function.\n\nIf that function is not needed, this one may return nothing to reduce memory allocations.\n\n\n\n\n\n","category":"function"},{"location":"implement/#SlidingDFTs.dataoffsets-implement","page":"Implementing SlidingDFTs","title":"SlidingDFTs.dataoffsets","text":"dataoffsets(method)\n\nReturn an integer or a vector of integers with the offsets of data samples that are needed by the given method to compute a sliding DFT.\n\nIf the code of SlidingDFTs.updatepdf! that dispatches on the type of method uses the function SlidingDFTs.previousdata, this function must return the integers that are used as the third argument (offset) of that function.\n\nIf that function is not needed (no past samples are used), this one may return nothing to reduce memory allocations.\n\n\n\n\n\n","category":"function"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"The fallback methods of those functions return nothing, as if neither previousdata or previousdft had to be used.","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"The implementation of updatepdf! given in the previous example does use previousdata - with the default offset value, equal to zero - so the following is required in this case:","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"SlidingDFTs.dataoffsets(::MyBasicSDFT) = 0","category":"page"},{"location":"implement/","page":"Implementing SlidingDFTs","title":"Implementing SlidingDFTs","text":"On the other hand there is no need to define SlidingDFTs.dftback in this case, since as it has been noted above, the interface of of updatedft! assumes that the most recent DFT is already contained in its first argument, so it is not necessary to use the function previousdft to get it.","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SlidingDFTs","category":"page"},{"location":"#SlidingDFTs","page":"Home","title":"SlidingDFTs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"SlidingDFTs is a Julia package to compute Sliding Discrete Fourer Transforms recursively, over one-dimensional series of values.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package is not yet registered. Add it with ]add https://github.com/heliosdrm/SlidingDFTs.jl in the \"pkg mode\" of the REPL, or in the \"standard REPL\":","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"https://github.com/heliosdrm/SlidingDFTs.jl\")","category":"page"},{"location":"#Basic-usage","page":"Home","title":"Basic usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package must be used with some implementation of AbstractFFTs such as FFTW, FastTransforms or RustFFT.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The basic Sliding Discrete Fourier Transform (SDFT) of a one-dimensional series of values x, using a window of length n, is calculated as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Step 1: Setup the method for a SDFT of length n:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using SlidingDFTs\nusing FFTW # or another package implementing AbstractFFTs\n\nmethod = SDFT(n)","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the API for other methods to compute sliding DFTs besides the basic SDFT.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Step 2: Create an iterator of the SDFT over x, with the function sdft. This is typically used in a loop:","category":"page"},{"location":"","page":"Home","title":"Home","text":"for spectrum in sdft(method, x)\n    # `spectrum` is a `Vector{Complex(eltype(x))}` of length `n`\nend","category":"page"},{"location":"#Improve-performance-with-\"unsafe\"-iterations","page":"Home","title":"Improve performance with \"unsafe\" iterations","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The iterator created by sdft produces a new vector of complex values at each iteration. For better performance, it is possible to call this function with the keyword argument safe set to false:","category":"page"},{"location":"","page":"Home","title":"Home","text":"iterator = sdft(method, x, safe=false)","category":"page"},{"location":"","page":"Home","title":"Home","text":"This iterator will produce only one vector that will be mutated at each iteration. Use this with caution, as any modification of the resulting vector may lead to unexpected results in subsequent iterations.","category":"page"},{"location":"#Using-SlidingDFTs-with-stateful-iterators","page":"Home","title":"Using SlidingDFTs with stateful iterators","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"By default, this package computes sliding DFTs traversing sequentially the data series x, which can be any kind of iterator. In the case of stateful iterators (i.e. those that are modified upon each iteration, like Base.Channels),  sdft(method, x) will also be a stateful iterator that will \"consume\" as many items of x as the length of the computed DFT in the first iteration, and one additional item in every subsequent iteration.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Apart from that consideration, it is safe to apply sliding DFTs to stateful iterators, since past samples of x already used in previous iterations, which are often required for the computations, are temporarily stored in an array — in internal variables that users do not need to deal with.","category":"page"}]
}