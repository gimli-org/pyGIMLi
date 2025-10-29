#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import shutil
import sys

try:
    import pygimli as pg
except:
    class pg:
        def _g(*args):
            print(*args)

        def warn(self, *args):
            print(*args)

        def error(self, *args):
            print(*args)

from environment_for_pygimli_build import settings

import hand_made_wrappers

from optparse import OptionParser
optionParser = OptionParser("usage: %prog [options]")
optionParser.add_option("", "--extra-includes", dest="extraIncludes")
optionParser.add_option("", "--extra-path", dest="extraPath")
optionParser.add_option("", "--caster", dest="caster")
optionParser.add_option("", "--clang", dest="clang", default="")

(options, args) = optionParser.parse_args()

print('##############')
print(options)
print(args)
print('##############')

if options.caster:
    settings.caster_path = options.caster

if options.extraPath:
    print("insert extra path: ", options.extraPath)
    sys.path.insert(0, options.extraPath)

from pygccxml import parser

import logging

from pygccxml import utils

# patch for >gcc-14
import patch_pygccxml

logger = utils.loggers.cxx_parser
logger.setLevel(logging.INFO)

from pygccxml import declarations
from pygccxml.declarations import access_type_matcher_t

from pyplusplus import code_creators, module_builder, messages, decl_wrappers
from pyplusplus.module_builder import call_policies
from pyplusplus.decl_wrappers.doc_extractor import doc_extractor_i

import hashlib

MAIN_NAMESPACE = 'GIMLI'


def samefile(sourcefile, destfile):
    """
    """
    if not os.path.exists(destfile):
        return False
    if not os.path.exists(sourcefile):
        return False

    inhash = hashlib.md5()
    inhash.update(open(sourcefile).read().encode('utf-8'))

    outhash = hashlib.md5()
    outhash.update(open(destfile).read().encode('utf-8'))

    if inhash.digest() != outhash.digest():
        return False

    # probably don't need these as hash should handle it fine..
    if os.stat(sourcefile).st_mtime > os.stat(destfile).st_mtime:
        return False
    if os.stat(sourcefile).st_size != os.stat(destfile).st_size:
        return False

    return True


class decl_starts_with (object):

    def __init__(self, prefix):
        self.prefix = prefix

    def __call__(self, decl):
        return self.prefix in decl.name


def exclude(method, return_type='', name='', symbol=''):
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
                # print("exclude name", f)
                f.exclude()

    for funct in symbol:
        if len(funct):
            fun = method(symbol=funct, allow_empty=True)

            for f in fun:
                # print("exclude symbol", f)
                f.exclude()


def setMemberFunctionCallPolicieByReturn(mb, MemberRetRef, callPolicie, verbose=False):
    for ref in MemberRetRef:
        memFuns = mb.global_ns.member_functions(
            return_type=ref,
            allow_empty=True)

        if verbose is True:
            print(ref, len(memFuns))

        for memFun in memFuns:
            memFun.call_policies = \
                call_policies.return_value_policy(callPolicie)
            # if memFun.call_policies:
            #     continue
            # else:


class docExtractor(doc_extractor_i):

    def __init__(self):
        doc_extractor_i.__init__(self)

    # def __call__(self, decl):
        # print "__call__(self, decl):", decl
        # print decl.location.file_name
        # print decl.location.line
        # return ""Doku here""

    def escape_doc(self, doc):
        return '"' + doc + '"'

    def extract(self, decl):
        # print(decl.location.file_name)
        # print(decl.location.line)
        # print("extract(self, decl):", decl)
        return "Doku coming soon"


def generate(defined_symbols, extraIncludes):
    messages.disable(
        messages.W1005  # using a non public variable type for arguments or returns
        , messages.W1006  # `Py++` need your
                          # help to expose function that takes > as argument/returns C++ arrays.
                          # Take a look on "Function Transformation" > functionality and define
                          # the transformation.
        , messages.W1007  # more than 10 args -> BOOST_PYTHON_MAX_ARITY is set
        , messages.W1009  # execution error W1009: The function takes as argument (name=pFunIdx, pos=1) >
                          # non-const reference to Python immutable type - function could not be called > from Python
        , messages.W1010  #: The function introduces registration order problem.
        , messages.W1014  # "operator*" is not supported. See
        , messages.W1016  # `Py++` does not exports non-const casting operators
        # Warnings 1020 - 1031 are all about why Py++ generates wrapper for class X
        , messages.W1023  # Py++` will generate class wrapper - there are few functions that should be
                          # redefined in class wrapper
        , messages.W1025  # `Py++` will generate class wrapper - class contains "c_" - T* > member variable
        , messages.W1026  # `Py++` will generate class wrapper - class contains "arr_" - T& > member variable
        , messages.W1027  # `Py++` will generate class wrapper - class contains "mat_" - > array member variable
        , messages.W1035  # error. `Py++` can not expose static pointer member variables.
        , messages.W1036  # error. `Py++` can not expose pointer to Python immutable > member variables. This
                          # could be changed in future.
        , messages.W1040  # error. The declaration is unexposed, but there are other > declarations, which
                          # refer to it. This could cause "no to_python converter > found" run
                          # time error
        , messages.W1042  #: `Py++` can not find out container value_type( mapped_type )
        # This is serious and lead to RuntimeError: `Py++` is going to write different content to the same file
        #, messages.W1047 # There are two or more classes that use same > alias("MatElement"). Duplicated aliases causes
                         # few problems, but the main one > is that some of the classes will not
                         # be exposed to Python.Other classes : >
        , messages.W1049  # This method could not be overriden in Python - method returns >
                          # reference to local variable!
        , messages.W1052  # `Py++` will not expose free operator
    )

    logger.debug("Install SRC:  " + os.path.abspath(__file__))
    logger.debug("Execute from: " + os.getcwd())

    sourcedir = os.path.dirname(os.path.abspath(__file__))
    sourceHeader = os.path.abspath(sourcedir + "/" + r"pygimli.h")
    gimliInclude = os.path.dirname(
                         os.path.abspath(sourcedir + "/../src/" + r"gimli.h"))
    settings.includesPaths.append(gimliInclude)

    xml_cached_fc = parser.create_cached_source_fc(
        sourceHeader, settings.module_name + '.cache')

    import platform

    defines = ['PYGIMLI_CAST', 'HAVE_BOOST_THREAD_HPP']
    caster = 'castxml'
    compiler_path = options.clang

    if platform.system() == 'Windows':
        if platform.architecture()[0] == '64bit':
            #compiler_path = 'C:/msys64/mingw64/bin/clang++'
            if sys.platform == 'darwin':
                pass
            else:
                defines.append('_WIN64')
                defines.append('MS_WIN64')

                logger.info('Marking win64 for gccxml')
        else:
            pass
            #compiler_path = 'C:/msys32/mingw32/bin/clang++'

    if len(compiler_path) == 0:
        compiler_path = None

    for define in [settings.gimli_defines, defined_symbols]:
        if len(define) > 0:
            defines.append(define)

    try:
        if sys.platform == 'win32':
            # os.name == 'nt' (default on my mingw)results in wrong commandline
            # for gccxml
            os.name = 'mingw'
            casterpath = settings.caster_path.replace('\\', '\\\\')
            casterpath = settings.caster_path.replace('/', '\\')

            if '.exe' not in casterpath:
                casterpath += '\\' + caster + '.exe'

            cflags = ''
        else:
            casterpath = settings.caster_path
            cflags = '-std=c++11'

    except Exception as e:
        logger.info("caster_path=%s" % casterpath)
        logger.info(str(e))
        raise Exception("Problems determine castxml binary")

    settings.includesPaths.insert(0, os.path.abspath(extraIncludes))

    logger.info("caster_path=%s" % casterpath)
    logger.info("working_directory=%s" % settings.gimli_path)
    logger.info("include_paths=%s" % settings.includesPaths)
    logger.info("define_symbols=%s" % defines)
    logger.info("compiler_path=%s" % compiler_path)
    logger.info("indexing_suite_version=2")

    xml_generator_config = parser.xml_generator_configuration_t(
                                        xml_generator=caster,
                                        xml_generator_path=casterpath,
                                        working_directory=settings.gimli_path,
                                        include_paths=settings.includesPaths,
                                        define_symbols=defines,
                                        ignore_gccxml_output=False,
                                        cflags=cflags,
                                        compiler_path=compiler_path)

    mb = module_builder.module_builder_t([xml_cached_fc],
                                         indexing_suite_version=2,
                                         xml_generator_config=xml_generator_config
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

        pg.error(e)
        sys.exit()

    logger.info("Apply custom rvalues.")

    # START manual r-value converters
    rvalue_converters = [
        'register_pysequence_to_indexarray_conversion',
        'register_pysequence_to_rvector_conversion',
        #'register_pysequence_to_cvector_conversion', .. will not work
        'register_pysequence_to_ivector_conversion',
        'register_pysequence_to_bvector_conversion',
        'register_pysequence_to_r3vector_conversion',
        'register_pysequence_to_StdVectorR3Vector_conversion',
        'register_pytuple_to_pos_conversion',
        'register_numpy_to_int64_conversion',
        'register_numpy_to_uint64_conversion',
        'register_numpy_to_double_conversion',
        'register_numpy_to_rmatrix_conversion',
        'register_numpy_to_rdensematrix_conversion',
        'register_numpy_to_int32_conversion',
        'register_numpy_to_uint32_conversion',
        'register_numpy_to_uint8_conversion',
        'register_numpy_to_uchar_conversion',
    ]
    for converter in rvalue_converters:
        mb.add_declaration_code('void %s();' % converter)
        mb.add_registration_code('%s();' % converter)
    # END manual r-value converters

    custom_rvalue_path = os.path.join(
        os.path.abspath(os.path.dirname(__file__)), 'custom_rvalue.cpp')

    logger.info("Starting to exclude stuff that we don't need "
                "or that is known to be spurious.")

    exclude(main_ns.variables,
            name=[
                'Triangle6_S1',
                'Triangle6_S2',
                'Triangle6_S3',
                'HexahedronFacesID',
                'Hexahedron20FacesID',
                'TetrahedronFacesID',
                'HexahedronSplit5TetID',
                'HexahedronSplit6TetID',
                'TriPrismFacesID',
                'TriPrimSplit3TetID',
                'NodeCoordinates',
                'EdgeCoordinates',
                'TriCoordinates',
                'QuadCoordinates',
                'TetCoordinates',
                'HexCoordinates',
                'PrismCoordinates',
                'PyramidCoordinates',
                'PyramidFacesID',
                'Tet10NodeSplit',
                'Tet10NodeSplitZienk',
                'Hex20NodeSplit',
                'Prism15NodeSplit',
                'Pyramid13NodeSplit',
                ]
            )

    exclude(main_ns.free_functions,
            return_type=[
                'float *',
                'float &',
                "::GIMLI::__VectorExpr< double, " +
                "GIMLI::__VectorUnaryExprOp< double, " +
                "GIMLI::VectorIterator< double >, GIMLI::ABS_ > >"],
            name=[
                'strReplaceBlankWithUnderscore',
                'toStr',
                'toInt',
                'toFloat',
                'toDouble',
                'str',
                'getRowSubstrings',
                'getNonEmptyRow',
                'getSubstrings',
                'abs',
                'type',
                '::GIMLI::print',
                'print',
                'range',
             ])

    exclude(main_ns.free_operators,
            name=[''],
            return_type=['::std::ostream &', '::std::istream &']
            )

    exclude(main_ns.classes,
            name=['ABS_', 'ACOT', 'ATAN', 'COS', 'COT', 'EXP', 'EXP10',
                  'ABS_', 'LOG', 'LOG10', 'SIGN', 'SIN', 'SQRT', 'SQR',
                  'TAN', 'TANH',
                  'PLUS', 'MINUS', 'MULT', 'DIVID', 'BINASSIGN', 'cerrPtr',
                  'ISINF', 'ISINFNAN',
                  'cerrPtrObject', 'coutPtr', 'coutPtrObject', 'deletePtr',
                  'edge_',
                  'distancePair_', 'IPCMessage', 'PythonGILSave',
                  'XAxis__', 'YAxis__', 'ZAxis__',
                'Variable',
                'Electrode',
                'ElectrodeShape',
                'ElectrodeShapeDomain',
                'ElectrodeShapeEntity',
                'ElectrodeShapeNode',
                'ElectrodeShapeNodesWithBypass',
                'FunctionDD',
                'H2SparseMapMatrix',
                #'::GIMLI::VectorIterator<GIMLI::Pos>', # needed
                #'::GIMLI::VectorIterator<bool>',
                #'::GIMLI::VectorIterator<double>',
                #'::GIMLI::VectorIterator<long>',
                # ::GIMLI::VectorIterator<unsigned long>', needed
                #'::GIMLI::VectorIterator<std::complex<double> >',
                  ]
            )

    exclude(main_ns.member_functions,
            name=['begin',
                  'end',
                  'val'],
            return_type=['']
            )

    exclude(main_ns.member_operators,
            symbol=[''])

    for f in main_ns.declarations:
        if isinstance(f, decl_wrappers.calldef_wrapper.free_function_t):
            if (str(f.return_type).find('GIMLI::VectorExpr') != -1):
                f.exclude()

    ex = ['::GIMLI::MatrixElement',
          '::GIMLI::__VectorUnaryExprOp',
          '::GIMLI::__VectorBinaryExprOp',
          '::GIMLI::__ValVectorExprOp',
          '::GIMLI::__VectorValExprOp',
          '::GIMLI::__VectorExpr',
          '::GIMLI::Expr',
          '::GIMLI::InversionBase',
          '::GIMLI::XAxis__',
          '::GIMLI::YAxis__',
          'GIMLI::ZAxis__',
          'GIMLI::ElectrodeShape',
          'GIMLI::MatrixElement',
          'GIMLI::__VectorUnaryExprOp',
          'GIMLI::__VectorBinaryExprOp',
          'GIMLI::__ValVectorExprOp',
          'GIMLI::__VectorValExprOp',
          'GIMLI::__VectorExpr',
          'GIMLI::Expr',
          'GIMLI::InversionBase',
          'std::vector<unsigned long',
          'std::vector<bool',
          'std::vector<double',
          'std::set<std',
        #   'std::set<Node*', # we need them
        #   'std::set<Cell*', # we need them
        #   'std::set<Boundary*', # we need them
          'Variable',
          'Variable<GIMLI::XAxis__',
          'FunctionDD',
          '::GIMLI::print',
          'print',
          'range',
          ]

    def excludeIfInDecl(c, excludes):
        for e in excludes:
            if c.decl_string.find(e) > -1:
                try:
                    c.exclude()
                    logger.debug("Exclude: " + str(c))
                    return True
                except BaseException as _:
                    logger.debug("Fail to exclude: " + str(c))

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
        except:
            pass

        try:
            for mem in c.variables():
                excludeIfInDecl(mem, ex)
        except:
            pass

        for mem in c.operators():
            excludeIfInDecl(mem, ex)

            ##########
            ### exclude all callable
            ##########
            if "operator()" in mem.name:

                #print("Exclude operator: " + str(mem))
                if not 'GIMLI::Trans' in str(mem):
                    logger.debug("Exclude callable for: " + str(mem))
                    mem.exclude()
                else:
                    pass
                    #print("skipping to exclude: " + str(mem))
                    #exit()

        # print('#'*100)
        # print(c, c.name)
        if c.name.startswith('Vector<unsigned long>'):
            # print('         ', c.name)
            for mem in c.constructors():
                # print("mem", mem, mem.decl_string)
                if mem.decl_string.find('( ::GIMLI::Index )') > -1:
                    logger.debug("Exclude: " + str(mem))
                    mem.exclude()


        if c.name.startswith("Matrix<double>") or \
            c.name.startswith("DenseMatrix<double>"):

            for mem in c.constructors():
                pass
                #print("\t constructors:", mem)
            for mem in c.operators():
                # print("\t operator:", mem)
                # print("\t\t\t:", mem.name)
                # print("\t\t\t:", mem.decl_string)

                if "operator()" in mem.name:
                    #print("Exclude: " + str(mem))
                    logger.debug("Exclude: " + str(mem))
                    mem.exclude()


    #sys.exit()

    try:
        mb.calldefs(access_type_matcher_t('protected')).exclude()
        mb.calldefs(access_type_matcher_t('private')).exclude()
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

    setMemberFunctionCallPolicieByReturn(mb,
        ['::std::string *', 'float *',
        'double *', 'double const *',
         'int *',
         'char *',
         'long *',
         'long int *',
         'long long int *',
         'unsigned int *',
         'long unsigned int *',
         'unsigned long long int *',
         'long long unsigned int *',
         '::GIMLI::Index *', '::GIMLI::SIndex *', 'bool *'],
        call_policies.return_pointee_value)

    setMemberFunctionCallPolicieByReturn(mb, ['::std::string &',
                                              'char &',
                                              'float &',
                                              'double &',
                                              'int &',
                                              'long &',
                                              'long int &',
                                              'long long int &',
                                              'unsigned int &',
                                              'long unsigned int &',
                                              'long long unsigned int &',
                                              'unsigned long long int &',
                                              '::GIMLI::Index &',
                                              '::GIMLI::SIndex &',
                                              'bool &'
                                              ], call_policies.return_by_value)

    setMemberFunctionCallPolicieByReturn(mb, ['::GIMLI::RVector &', #0
                                              '::GIMLI::RVector const &', #23
                                              '::GIMLI::CVector &', #0
                                              '::GIMLI::CVector const &', #0
                                              '::GIMLI::IVector &', #0
                                              '::GIMLI::IVector const &', #1
                                              '::GIMLI::IndexArray &', #0
                                              '::GIMLI::IndexArray const &', #1
                                              '::GIMLI::SVector &', #0
                                              '::GIMLI::SVector const &', #0
                                              '::GIMLI::RMatrix &', #5
                                              '::GIMLI::RMatrix const &', #5
                                              '::GIMLI::CMatrix &', #0
                                              '::GIMLI::CMatrix const &', #0
                                              '::GIMLI::RSmallMatrix &', #0
                                              '::GIMLI::RSmallMatrix const &', #1
                                              '::GIMLI::RDenseMatrix &', #0
                                              '::GIMLI::RDenseMatrix const &', #0
                                              ],
                                        call_policies.reference_existing_object, verbose=False)


    # setMemberFunctionCallPolicieByReturn(mb,
    #                                   ['::GIMLI::VectorIterator<double> &']
    # , call_policies.copy_const_reference)
    # setMemberFunctionCallPolicieByReturn(mb, [
    # ,  'double &' ]
    # , call_policies.reference_existing_object)

    # call_policies.return_value_policy(call_policies.reference_existing_object)
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
            if not mem_fun.call_policies and \
                    (declarations.is_reference(mem_fun.return_type) or
                        declarations.is_pointer(mem_fun.return_type)):
                # print mem_fun
                # mem_fun.exclude()
                mem_fun.call_policies = call_policies.return_value_policy(
                    call_policies.reference_existing_object)
                # mem_fun.call_policies = \
                #    call_policies.return_value_policy(call_policies.return_pointee_value)
                # mem_fun.call_policies = \
                #    call_policies.return_value_policy(call_policies.return_opaque_pointer)
                # mem_fun.call_policies = \
                #   call_policies.return_value_policy(call_policies.copy_non_const_reference)

    logger.info("Create API documentation from Doxygen comments.")
    # Now it is the time to give a name to our module
    from doxygen import doxygen_doc_extractor
    extractor = doxygen_doc_extractor()

    logger.info("Create code creator.")
    mb.build_code_creator(settings.module_name, doc_extractor=extractor)

    # It is common requirement in software world-each file should have license
    # mb.code_creator.license = '//Boost Software
    # License(http://boost.org/more/license_info.html)'

    # I don't want absolute includes within code
    mb.code_creator.user_defined_directories.append(os.path.abspath('.'))

    # And finally we can write code to the disk
    def ignore(val):
        pass
    logger.info("Create bindings code.")
    mb.split_module('./generated', on_unused_file_found=ignore)

    additional_files = [
        os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'custom_rvalue.cpp'),
        os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'explicit_instances.cpp'),
        os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'generators.h'),
        os.path.join(
            os.path.abspath(os.path.dirname(__file__)), 'tuples.hpp')]

    logger.info("Add additional files.")

    for sourcefile in additional_files:
        p, filename = os.path.split(sourcefile)
        destfile = os.path.join('./generated', filename)

        if not samefile(sourcefile, destfile):
            shutil.copy(sourcefile, './generated')
            logger.info("Updated " + filename +
                        "as it was missing or out of date")

if __name__ == '__main__':

    defined_symbols = ''

    generate(defined_symbols, options.extraIncludes)
