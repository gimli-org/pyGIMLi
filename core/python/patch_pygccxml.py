import pygccxml

from pygccxml.parser import project_reader_t

# for gcc-14 this is needed
def __patch_relink_declarated_types__(self, leaved_classes, declarated_types):
    mangled_leaved_classes = {}
    for leaved_class in leaved_classes.values():
        mangled_leaved_classes[
            self._create_name_key(leaved_class)] = leaved_class

    for decl_wrapper_type in declarated_types:
        # it is possible, that cache contains reference to dropped class
        # We need to clear it
        decl_wrapper_type.cache.reset()
        if isinstance(
                decl_wrapper_type.declaration,
                pygccxml.declarations.class_t):
            key = self._create_key(decl_wrapper_type.declaration)
            if key in leaved_classes:
                decl_wrapper_type.declaration = leaved_classes[key]
            else:

                name = decl_wrapper_type.declaration._name
                if name == "":
                    # Happens with gcc5, castxml + std=c++11
                    # See issue #45
                    continue
                if name.startswith("__vmi_class_type_info_pseudo"):
                    continue
                if name == "rebind<std::__tree_node" + \
                        "<std::basic_string<char>, void *> >":
                    continue

                ########## patch on
                print('+'*80)
                print('key',key)
                print('_'*80)
                print('name', name)
                print('-'*80)
                continue
                ########## patch off

                msg = []
                msg.append(
                    "Unable to find out actual class definition: '%s'." %
                    decl_wrapper_type.declaration._name)
                msg.append((
                    "Class definition has been changed from one " +
                    "compilation to an other."))
                msg.append((
                    "Why did it happen to me? Here is a short list " +
                    "of reasons: "))
                msg.append((
                    "    1. There are different preprocessor " +
                    "definitions applied on same file during compilation"))
                msg.append("    2. Bug in pygccxml.")
                raise Exception(os.linesep.join(msg))
        elif isinstance(
                decl_wrapper_type.declaration,
                pygccxml.declarations.class_declaration_t):
            key = self._create_name_key(decl_wrapper_type.declaration)
            if key in mangled_leaved_classes:
                decl_wrapper_type.declaration = mangled_leaved_classes[key]

project_reader_t._relink_declarated_types = __patch_relink_declarated_types__
