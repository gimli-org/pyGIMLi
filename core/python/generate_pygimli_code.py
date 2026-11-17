#!/usr/bin/env python
"""Generates the pygimli python bindings using py++ and pygccxml."""

import hashlib
import logging
import os
import shutil
import sys
import tempfile
from optparse import OptionParser
from pathlib import Path

import hand_made_wrappers
from environment_for_pygimli_build import Settings


class Pg:
    """Print helper for the messages of the code generator."""

    def _g(*args):
        print(*args)

    def warn(self, *args):
        """Print a warning message."""
        print(*args)

    def error(self, *args):
        """Print an error message."""
        print(*args)


optionParser = OptionParser("usage: %prog [options]")
optionParser.add_option("", "--extra-includes", dest="extraIncludes")
optionParser.add_option("", "--extra-path", dest="extraPath")
optionParser.add_option("", "--castxml", dest="caster")
optionParser.add_option("", "--clang", dest="clang", default="")

(options, args) = optionParser.parse_args()

print("##############")
print(options)
print(args)
print("##############")

if options.caster:
    Settings.caster_path = options.caster

if options.extraPath:
    print("insert extra path: ", options.extraPath)
    sys.path.insert(0, options.extraPath)

# these imports must stay below the --extra-path handling: the option may
# add the location of pygccxml/pyplusplus to sys.path
from pygccxml import parser  # noqa: E402
from pygccxml import utils  # noqa: E402

# patch for >gcc-14 (applied by importing the module)
import patch_pygccxml  # noqa: E402,F401

from pygccxml import declarations  # noqa: E402
from pygccxml.declarations import (  # noqa: E402
    access_type_matcher_t,
)
from pyplusplus import module_builder, messages, decl_wrappers  # noqa: E402
from pyplusplus.module_builder import call_policies  # noqa: E402
from pyplusplus.decl_wrappers.doc_extractor import (  # noqa: E402
    doc_extractor_i,
)

# patch for >gcc-14
logger = utils.loggers.cxx_parser
logger.setLevel(logging.INFO)

MAIN_NAMESPACE = "GIMLI"


class CastXmlSetupError(RuntimeError):
    """Raised when the castxml binary can not be determined."""


def samefile(sourcefile, destfile):
    """Return True if destfile has the same content as sourcefile."""
    if not Path(destfile).exists() or not Path(sourcefile).exists():
        return False

    inhash = hashlib.md5(usedforsecurity=False)
    inhash.update(Path(sourcefile).read_bytes())

    outhash = hashlib.md5(usedforsecurity=False)
    outhash.update(Path(destfile).read_bytes())

    if inhash.digest() != outhash.digest():
        return False

    sourceStat = Path(sourcefile).stat()
    destStat = Path(destfile).stat()

    # probably don't need these as hash should handle it fine..
    return (
        sourceStat.st_mtime <= destStat.st_mtime
        and sourceStat.st_size == destStat.st_size
    )


class DeclStartsWith:
    """Filter declarations whose name contains a given prefix."""

    def __init__(self, prefix):
        self.prefix = prefix

    def __call__(self, decl):
        """Return True if the name of decl contains the prefix."""
        return self.prefix in decl.name


def exclude(method, return_type="", name="", symbol=""):
    """Exclude matching declarations from the generated bindings."""
    for funct in return_type:
        if len(funct):
            fun = method(return_type=funct, allow_empty=True)

            for f in fun:
                # print("exclude return type", f)
                f.exclude()

    for funct in name:
        if len(funct):
            fun = method(name=funct, allow_empty=True)

            for f in fun:
                # print("exclude "_", f)
                f.exclude()

    for funct in symbol:
        if len(funct):
            fun = method(symbol=funct, allow_empty=True)

            for f in fun:
                # print("exclude symbol", f)
                f.exclude()


def setMemberFunctionCallPolicieByReturn(
    mb, MemberRetRef, callPolicie, verbose=False
):
    """Set the return value policy for member functions by return type."""
    for ref in MemberRetRef:
        memFuns = mb.global_ns.member_functions(
            return_type=ref, allow_empty=True
        )

        if verbose is True:
            print(ref, len(memFuns))

        for memFun in memFuns:
            memFun.call_policies = call_policies.return_value_policy(
                callPolicie
            )
            # if memFun.call_policies:
            #     continue
            # else:


class DocExtractor(doc_extractor_i):
    """Doxygen doc extractor returning a fixed placeholder text."""

    def __init__(self):
        doc_extractor_i.__init__(self)

    # def __call__(self, decl):
    # print "__call__(self, decl):", decl
    # print decl.location.file_name
    # print decl.location.line
    # return ""Doku here""

    def escape_doc(self, doc):
        """Escape the doc string by wrapping it in double quotes."""
        return '"' + doc + '"'

    def extract(self, decl):
        """Return the documentation for a declaration."""
        # print(decl.location.file_name)
        # print(decl.location.line)
        # print("extract(self, decl):", decl)
        return "Doku coming soon"


def generate(defined_symbols, extraIncludes):
    """Generate the python binding sources into ./generated."""
    messages.disable(
        messages.W1005,  # using a non public variable type for args or returns
        messages.W1006,  # `Py++` needs help to expose functions that
        # take/return C++ arrays (see "Function
        # Transformation" and define the transformation)
        messages.W1007,  # more than 10 args -> BOOST_PYTHON_MAX_ARITY
        messages.W1009,  # the function takes a non-const reference to a
        # Python immutable type as argument (pFunIdx) ->
        # it could not be called from Python
        messages.W1010,  #: The function introduces registration order
        # problem.
        messages.W1014,  # "operator*" is not supported. See
        messages.W1016,  # `Py++` does not exports non-const casting operators
        # Warnings 1020 - 1031 are all about why Py++ generates wrappers
        # for class X
        messages.W1023,  # `Py++` will generate a class wrapper - there are
        # few functions that should be redefined in it
        messages.W1025,  # `Py++` will generate a class wrapper - the class
        # contains a "c_" - T* member variable
        messages.W1026,  # `Py++` will generate a class wrapper - the class
        # contains an "arr_" - T& member variable
        messages.W1027,  # `Py++` will generate a class wrapper - the class
        # contains a "mat_" - array member variable
        messages.W1035,  # `Py++` can not expose static pointer members
        messages.W1036,  # `Py++` can not expose pointers to Python
        # immutable member variables (could change in future)
        messages.W1040,  # the declaration is unexposed, but other
        # declarations refer to it, which could cause a
        # "no to_python converter found" runtime error
        messages.W1042,  #: `Py++` can not find container value_type
        # This is serious and leads to RuntimeError: `Py++` is going to
        # write different content to the same file
        # , messages.W1047  # two or more classes use the same alias
        #                    # ("MatElement"). Duplicated aliases cause few
        #                    # problems, but some classes will then not be
        #                    # exposed to Python
        messages.W1049,  # this method could not be overridden in Python -
        # it returns a reference to a local variable
        messages.W1052,  # `Py++` will not expose free operator
    )

    logger.debug("Install SRC:  " + str(Path(__file__).resolve()))
    logger.debug("Execute from: " + str(Path.cwd()))

    sourceDir = Path(__file__).resolve().parent
    sourceHeader = str((sourceDir / r"pygimli.h").resolve())
    gimliInclude = str((sourceDir / r"../src/gimli.h").resolve().parent)
    Settings.includes_paths.append(gimliInclude)

    xml_cached_fc = parser.create_cached_source_fc(
        sourceHeader, Settings.module_name + ".cache"
    )

    import platform

    defines = ["PYGIMLI_CAST", "HAVE_BOOST_THREAD_HPP"]
    caster = "castxml"
    compiler_path = options.clang

    if platform.system() == "Windows":
        if platform.architecture()[0] == "64bit":
            # compiler_path = 'C:/msys64/mingw64/bin/clang++'
            if sys.platform == "darwin":
                pass
            else:
                defines.append("_WIN64")
                defines.append("MS_WIN64")

                logger.info("Marking win64 for gccxml")
        else:
            pass
            # compiler_path = 'C:/msys32/mingw32/bin/clang++'

    if len(compiler_path) == 0:
        compiler_path = None

    for define in [Settings.gimli_defines, defined_symbols]:
        if len(define) > 0:
            defines.append(define)

    try:
        if sys.platform == "win32":
            # os.name == 'nt' (default on my mingw)results in wrong commandline
            # for gccxml
            os.name = "mingw"
            casterpath = Settings.caster_path.replace("\\", "\\\\")
            casterpath = Settings.caster_path.replace("/", "\\")

            if ".exe" not in casterpath:
                casterpath += "\\" + caster + ".exe"
            # GCC 16 libstdc++ defines __int128 specializations *both* via
            # __GLIBCXX_TYPE_INT_N_0 (the GCC built-in alias) *and* directly
            # as __int128.  When CastXML/Clang parses with both macros active
            # the two sets collide.  Undefining the GCC alias macro leaves only
            # the direct __int128 specialisations and resolves all
            # "redefinition of '__is_integral_helper<__int128>'" errors.
            cflags = (
                "-std=gnu++20 -U__GLIBCXX_TYPE_INT_N_0 "
                "-U__GLIBCXX_BITSIZE_INT_N_0"
            )
        else:
            casterpath = Settings.caster_path
            cflags = "-std=c++11"

    except Exception as e:
        logger.info(f"caster_path={casterpath}")
        logger.info(str(e))
        raise CastXmlSetupError("Problems determine castxml binary") from e

    Settings.includes_paths.insert(0, str(Path(extraIncludes).resolve()))

    # Workaround for castxml (clang 18) + GCC 16 headers incompatibility.
    # GCC 16 added [[__gnu__::__always_inline__]] to the friend template
    # operator- in __normal_iterator with a trailing decltype return type.
    # Clang 18 evaluates that decltype before the class body is complete,
    # producing "member access into incomplete type" errors for every
    # std::vector/std::string use deep inside the STL headers. Fix by
    # placing a patched copy of bits/stl_iterator.h (with the offending
    # attribute removed) first in the include search path.
    if sys.platform != "win32":
        _candidates = sorted(
            str(p)
            for p in Path("/usr/lib/gcc").glob(
                "*/*/include/g++*/bits/stl_iterator.h"
            )
        )
        # Only apply when GCC 16 headers are present (the __always_inline__
        # friend operator- regression was introduced in GCC 16).
        _candidates = [
            p
            for p in _candidates
            if "/gcc/x86_64" in p.replace("\\", "/")
            and int(p.split("/gcc/")[1].split("/")[1].split(".")[0]) >= 16
        ]
        if _candidates:
            _stl_iter_orig = _candidates[-1]  # newest GCC >= 16 version
            _patch_dir = tempfile.mkdtemp(prefix="gimli_castxml_compat_")
            _bits_dir = Path(_patch_dir) / "bits"
            _bits_dir.mkdir()
            with Path(_stl_iter_orig).open() as _f:
                _content = _f.read()
            # Replace the decltype trailing return type with explicit
            # difference_type to avoid clang 18 evaluating __lhs.base()
            # while __normal_iterator is still incomplete (GCC 16 regression).
            import re

            _content = re.sub(
                r"\[\[__nodiscard__,\s*__gnu__::__always_inline__\]\]"
                r"(\s+friend\s+// _GLIBCXX_RESOLVE_LIB_DEFECTS"
                r"\s+// 685[^\n]+\s+constexpr\s+)auto"
                r"(\s+operator-\([^)]+\)\s+noexcept)"
                r"\s*\n\s*->\s*decltype\(__lhs\.base\(\)\s*-\s*__rhs\.base\(\)\)",
                r"[[__nodiscard__]]\1difference_type\2",
                _content,
            )
            with (_bits_dir / "stl_iterator.h").open("w") as _f:
                _f.write(_content)
            Settings.includes_paths.insert(0, _patch_dir)
            logger.info(
                "castxml compat: patched %s -> %s",
                _stl_iter_orig,
                str(_bits_dir / "stl_iterator.h"),
            )

    logger.info(f"caster_path={casterpath}")
    logger.info(f"working_directory={Settings.gimli_path}")
    logger.info(f"include_paths={Settings.includes_paths}")
    logger.info(f"define_symbols={defines}")
    logger.info(f"compiler_path={compiler_path}")
    logger.info("indexing_suite_version=2")

    xml_generator_config = parser.xml_generator_configuration_t(
        xml_generator=caster,
        xml_generator_path=casterpath,
        working_directory=Settings.gimli_path,
        include_paths=Settings.includes_paths,
        define_symbols=defines,
        ignore_gccxml_output=False,
        cflags=cflags,
        compiler_path=compiler_path,
    )

    mb = module_builder.module_builder_t(
        [xml_cached_fc],
        indexing_suite_version=2,
        xml_generator_config=xml_generator_config,
    )

    logger.info("Reading of c++ sources done.")

    mb.classes().always_expose_using_scope = True
    mb.calldefs().create_with_signature = True

    global_ns = mb.global_ns
    global_ns.exclude()
    main_ns = global_ns.namespace(MAIN_NAMESPACE)
    main_ns.include()

    # for c in main_ns.free_functions():
    #   print(c)
    #   if 'pow' in c.name:
    #       print(c)
    #       print(c.name)
    #   if c.decl_string.startswith('::GIMLI::pow'):
    #       print(c)
    #       print(c.name)
    # sys.exit()

    logger.info("Apply handmade wrappers.")
    try:
        hand_made_wrappers.apply(mb)
    except BaseException as e:
        import traceback

        traceback.print_exc(file=sys.stdout)

        Pg.error(e)
        sys.exit()

    logger.info("Apply custom rvalues.")

    # START manual r-value converters
    rvalue_converters = [
        "register_pysequence_to_indexarray_conversion",
        "register_pysequence_to_rvector_conversion",
        #'register_pysequence_to_cvector_conversion', .. will not work
        "register_pysequence_to_ivector_conversion",
        "register_pysequence_to_bvector_conversion",
        "register_pysequence_to_r3vector_conversion",
        "register_pysequence_to_StdVectorR3Vector_conversion",
        "register_pytuple_to_pos_conversion",
        "register_numpy_to_int64_conversion",
        "register_numpy_to_uint64_conversion",
        "register_numpy_to_double_conversion",
        "register_int_to_double_conversion",
        "register_numpy_to_rmatrix_conversion",
        "register_numpy_to_rdensematrix_conversion",
        "register_numpy_to_int32_conversion",
        "register_numpy_to_uint32_conversion",
        "register_numpy_to_uint8_conversion",
        "register_numpy_to_uchar_conversion",
    ]
    for converter in rvalue_converters:
        mb.add_declaration_code(f"void {converter}();")
        mb.add_registration_code(f"{converter}();")
    # END manual r-value converters

    logger.info(
        "Starting to exclude stuff that we don't need "
        "or that is known to be spurious."
    )

    exclude(
        main_ns.variables,
        name=[
            "Triangle6_S1",
            "Triangle6_S2",
            "Triangle6_S3",
            "HexahedronFacesID",
            "Hexahedron20FacesID",
            "TetrahedronFacesID",
            "HexahedronSplit5TetID",
            "HexahedronSplit6TetID",
            "TriPrismFacesID",
            "TriPrimSplit3TetID",
            "NodeCoordinates",
            "EdgeCoordinates",
            "TriCoordinates",
            "QuadCoordinates",
            "TetCoordinates",
            "HexCoordinates",
            "PrismCoordinates",
            "PyramidCoordinates",
            "PyramidFacesID",
            "Tet10NodeSplit",
            "Tet10NodeSplitZienk",
            "Hex20NodeSplit",
            "Prism15NodeSplit",
            "Pyramid13NodeSplit",
        ],
    )

    exclude(
        main_ns.free_functions,
        return_type=[
            "float *",
            "float &",
            "::GIMLI::__VectorExpr< double, "
            + "GIMLI::__VectorUnaryExprOp< double, "
            + "GIMLI::VectorIterator< double >, GIMLI::ABS_ > >",
        ],
        name=[
            "strReplaceBlankWithUnderscore",
            "toStr",
            "toInt",
            "toFloat",
            "toDouble",
            "str",
            "getRowSubstrings",
            "getNonEmptyRow",
            "getSubstrings",
            "abs",
            "type",
            "::GIMLI::print",
            "print",
            "range",
        ],
    )

    exclude(
        main_ns.free_operators,
        name=[""],
        return_type=["::std::ostream &", "::std::istream &"],
    )

    exclude(
        main_ns.classes,
        name=[
            "ABS_",
            "ACOT",
            "ATAN",
            "COS",
            "COT",
            "EXP",
            "EXP10",
            "ABS_",
            "LOG",
            "LOG10",
            "SIGN",
            "SIN",
            "SQRT",
            "SQR",
            "TAN",
            "TANH",
            "PLUS",
            "MINUS",
            "MULT",
            "DIVID",
            "BINASSIGN",
            "cerrPtr",
            "ISINF",
            "ISINFNAN",
            "cerrPtrObject",
            "coutPtr",
            "coutPtrObject",
            "deletePtr",
            "edge_",
            "distancePair_",
            "IPCMessage",
            "PythonGILSave",
            "XAxis__",
            "YAxis__",
            "ZAxis__",
            "Variable",
            "Electrode",
            "ElectrodeShape",
            "ElectrodeShapeDomain",
            "ElectrodeShapeEntity",
            "ElectrodeShapeNode",
            "ElectrodeShapeNodesWithBypass",
            "FunctionDD",
            "H2SparseMapMatrix",
            #'::GIMLI::VectorIterator<GIMLI::Pos>', # needed
            #'::GIMLI::VectorIterator<bool>',
            #'::GIMLI::VectorIterator<double>',
            #'::GIMLI::VectorIterator<long>',
            # ::GIMLI::VectorIterator<unsigned long>', needed
            #'::GIMLI::VectorIterator<std::complex<double> >',
        ],
    )

    exclude(
        main_ns.member_functions, name=["begin", "end", "val"], return_type=[""]
    )

    exclude(main_ns.member_operators, symbol=[""])

    for f in main_ns.declarations:
        if (
            isinstance(f, decl_wrappers.calldef_wrapper.free_function_t)
            and str(f.return_type).find("GIMLI::VectorExpr") != -1
        ):
            f.exclude()

    ex = [
        "::GIMLI::MatrixElement",
        "::GIMLI::__VectorUnaryExprOp",
        "::GIMLI::__VectorBinaryExprOp",
        "::GIMLI::__ValVectorExprOp",
        "::GIMLI::__VectorValExprOp",
        "::GIMLI::__VectorExpr",
        "::GIMLI::Expr",
        "::GIMLI::InversionBase",
        "::GIMLI::XAxis__",
        "::GIMLI::YAxis__",
        "GIMLI::ZAxis__",
        "GIMLI::ElectrodeShape",
        "GIMLI::MatrixElement",
        "GIMLI::__VectorUnaryExprOp",
        "GIMLI::__VectorBinaryExprOp",
        "GIMLI::__ValVectorExprOp",
        "GIMLI::__VectorValExprOp",
        "GIMLI::__VectorExpr",
        "GIMLI::Expr",
        "GIMLI::InversionBase",
        "std::vector<unsigned long",
        "std::vector<bool",
        "std::vector<double",
        "std::set<std",
        #   'std::set<Node*', # we need them
        #   'std::set<Cell*', # we need them
        #   'std::set<Boundary*', # we need them
        "Variable",
        "Variable<GIMLI::XAxis__",
        "FunctionDD",
        "::GIMLI::print",
        "print",
        "range",
    ]

    def excludeIfInDecl(c, excludes):
        for e in excludes:
            if c.decl_string.find(e) > -1:
                try:
                    c.exclude()
                except BaseException as _:
                    logger.debug("Fail to exclude: " + str(c))
                else:
                    logger.debug("Exclude: " + str(c))
                    return True

        return False

    for c in main_ns.free_functions():
        excludeIfInDecl(c, ex)

    for c in main_ns.classes():
        if excludeIfInDecl(c, ex) is True:
            continue

        for mem in c.constructors():
            excludeIfInDecl(mem, ex)

        try:
            for mem in c.member_functions():
                excludeIfInDecl(mem, ex)
        except Exception:
            pass

        try:
            for mem in c.variables():
                excludeIfInDecl(mem, ex)
        except Exception:
            pass

        for mem in c.operators():
            excludeIfInDecl(mem, ex)

            ##########
            ### exclude all callable
            ##########
            if "operator()" in mem.name:
                # print("Exclude operator: " + str(mem))
                if "GIMLI::Trans" not in str(mem):
                    logger.debug("Exclude callable for: " + str(mem))
                    mem.exclude()
                else:
                    pass
                    # print("skipping to exclude: " + str(mem))
                    # exit()

        # print('#'*100)
        # print(c, c.name)
        if c.name.startswith("Vector<unsigned long>"):
            # print('         ', c.name)
            for mem in c.constructors():
                # print("mem", mem, mem.decl_string)
                if mem.decl_string.find("( ::GIMLI::Index )") > -1:
                    logger.debug("Exclude: " + str(mem))
                    mem.exclude()

        if c.name.startswith("Matrix<double>") or c.name.startswith(
            "DenseMatrix<double>"
        ):
            for _mem in c.constructors():
                pass
                # print("\t constructors:", _mem)
            for mem in c.operators():
                # print("\t operator:", mem)
                # print("\t\t\t:", mem.name)
                # print("\t\t\t:", mem.decl_string)

                if "operator()" in mem.name:
                    # print("Exclude: " + str(mem))
                    logger.debug("Exclude: " + str(mem))
                    mem.exclude()

    # sys.exit()

    try:
        mb.calldefs(access_type_matcher_t("protected")).exclude()
        mb.calldefs(access_type_matcher_t("private")).exclude()
    except BaseException as _:
        pass

    # setMemberFunctionCallPolicieByReturn(mb, [ '::GIMLI::Node &'
    # , '::GIMLI::Cell &'
    # , '::GIMLI::Boundary &'
    # , '::GIMLI::Shape &'
    # , '::GIMLI::Node *'
    # , '::GIMLI::Cell *'
    # , '::GIMLI::Boundary *'
    # , '::GIMLI::Shape *'
    # ]
    # , call_policies.reference_existing_object)

    setMemberFunctionCallPolicieByReturn(
        mb,
        [
            "::std::string *",
            "float *",
            "double *",
            "double const *",
            "int *",
            "char *",
            "long *",
            "long int *",
            "long long int *",
            "unsigned int *",
            "long unsigned int *",
            "unsigned long long int *",
            "long long unsigned int *",
            "::GIMLI::Index *",
            "::GIMLI::SIndex *",
            "bool *",
        ],
        call_policies.return_pointee_value,
    )

    setMemberFunctionCallPolicieByReturn(
        mb,
        [
            "::std::string &",
            "char &",
            "float &",
            "double &",
            "int &",
            "long &",
            "long int &",
            "long long int &",
            "unsigned int &",
            "long unsigned int &",
            "long long unsigned int &",
            "unsigned long long int &",
            "::GIMLI::Index &",
            "::GIMLI::SIndex &",
            "bool &",
        ],
        call_policies.return_by_value,
    )

    setMemberFunctionCallPolicieByReturn(
        mb,
        [
            "::GIMLI::RVector &",  # 0
            "::GIMLI::RVector const &",  # 23
            "::GIMLI::CVector &",  # 0
            "::GIMLI::CVector const &",  # 0
            "::GIMLI::IVector &",  # 0
            "::GIMLI::IVector const &",  # 1
            "::GIMLI::IndexArray &",  # 0
            "::GIMLI::IndexArray const &",  # 1
            "::GIMLI::SVector &",  # 0
            "::GIMLI::SVector const &",  # 0
            "::GIMLI::RMatrix &",  # 5
            "::GIMLI::RMatrix const &",  # 5
            "::GIMLI::CMatrix &",  # 0
            "::GIMLI::CMatrix const &",  # 0
            "::GIMLI::RSmallMatrix &",  # 0
            "::GIMLI::RSmallMatrix const &",  # 1
            "::GIMLI::RDenseMatrix &",  # 0
            "::GIMLI::RDenseMatrix const &",  # 0
        ],
        call_policies.reference_existing_object,
        verbose=False,
    )

    # setMemberFunctionCallPolicieByReturn(mb,
    #                                   ['::GIMLI::VectorIterator<double> &']
    # , call_policies.copy_const_reference)
    # setMemberFunctionCallPolicieByReturn(mb, [
    # ,  'double &' ]
    # , call_policies.reference_existing_object)

    # call_policies.return_value_policy(
    #     call_policies.reference_existing_object)
    # call_policies.return_value_policy(call_policies.copy_non_const_reference)
    # call_policies.return_value_policy(call_policies.copy_const_reference)

    # addAutoConversions(mb)

    # excludeMemberByReturn(main_ns, ['::DCFEMLib::SparseMatrix<double> &'])
    # fun = mb.global_ns.member_functions('begin', allow_empty=True)
    # for f in fun:
    # f.exclude()

    # excludeFreeFunctionsByName(main_ns, ['strReplaceBlankWithUnderscore'
    # 'toStr', 'toInt', 'toFloat', 'toDouble',
    # 'getRowSubstrings', 'getNonEmptyRow', 'getSubstrings' ])

    # excludeFreeFunctionsByReturn(main_ns, [ 'float *', 'float &' ])
    # fun = ns.free_operators(return_type=funct, allow_empty=True)

    # excludeMemberOperators(main_ns, ['++', '--', '*'])

    # exclude all that does not match any predefined callpolicie

    excludeRest = True

    if excludeRest:
        mem_funs = mb.calldefs()

        for mem_fun in mem_funs:
            if mem_fun.call_policies:
                continue
            if not mem_fun.call_policies and (
                declarations.is_reference(mem_fun.return_type)
                or declarations.is_pointer(mem_fun.return_type)
            ):
                # print mem_fun
                # mem_fun.exclude()
                mem_fun.call_policies = call_policies.return_value_policy(
                    call_policies.reference_existing_object
                )
                # mem_fun.call_policies = \
                #    call_policies.return_value_policy(
                #        call_policies.return_pointee_value)
                # mem_fun.call_policies = \
                #    call_policies.return_value_policy(
                #        call_policies.return_opaque_pointer)
                # mem_fun.call_policies = \
                #   call_policies.return_value_policy(
                #       call_policies.copy_non_const_reference)

    logger.info("Create API documentation from Doxygen comments.")
    # Now it is the time to give a name to our module
    from doxygen import DoxygenDocExtractor

    extractor = DoxygenDocExtractor()

    logger.info("Create code creator.")
    mb.build_code_creator(Settings.module_name, doc_extractor=extractor)

    # It is common requirement in software world-each file should have license
    # mb.code_creator.license = '//Boost Software
    # License(http://boost.org/more/license_info.html)'

    # I don't want absolute includes within code
    mb.code_creator.user_defined_directories.append(str(Path.cwd()))

    # And finally we can write code to the disk
    def ignore(val):
        pass

    logger.info("Create bindings code.")
    mb.split_module("./generated", on_unused_file_found=ignore)

    moduleDir = Path(__file__).resolve().parent
    additional_files = [
        moduleDir / "custom_rvalue.cpp",
        moduleDir / "explicit_instances.cpp",
        moduleDir / "generators.h",
        moduleDir / "tuples.hpp",
    ]

    logger.info("Add additional files.")

    for sourcefile in additional_files:
        filename = sourcefile.name
        destfile = Path("./generated") / filename

        if not samefile(sourcefile, destfile):
            shutil.copy(sourcefile, "./generated")
            logger.info(
                "Updated " + filename + "as it was missing or out of date"
            )


if __name__ == "__main__":
    defined_symbols = ""

    generate(defined_symbols, options.extraIncludes)
