"""
Sphinx extension to add ReadTheDocs-style "Edit on GitHub" links to the
sidebar.

Loosely based on https://github.com/astropy/astropy/pull/347

| https://gist.github.com/westurner/d89c1ea1af05c5c514f9
| Src: https://gist.github.com/d89c1ea1af05c5c514f9.git

"""

import collections
import itertools
import json
import logging
import os
import re
import sys
try:
    import urlparse
except:
    import urllib.parse as urlparse
import warnings


__licence__ = 'BSD (3 clause)'

log = logging.getLogger(__name__)


class Context(object):
    type_ = None
    src_type = None         # 'git', 'hg'
    src_https_tmpl = None   # 'https://github.com/{project}'
    src_url_tmpl = None     #  src_https_tmpl
    src_ssh_tmpl = None     # 'ssh://git@github.com/{project}'
    src_native_tmpl = None  # 'git://github.com/{project}'
    src_view_url_tmpl = '{project_url}/{view}/{rev}/{src_path}{path}'
    src_clonecmd_tmpl = '{src_type} clone'

    src_issue_regexes = [
        r"#[\d]+",  # TODO: {,10}
    ]

    def __init__(
            self,
            app=None,
            pagename=None,
            templatename=None,
            context=None,
            doctree=None,
            path=None,
            project=None):
        """
        Create a new Context from
        urlstr and/or more specifically, project
        """
        if app:
            project = app.config.srclink_project
        if project is None:
            raise ValueError("project is None")
        self.project = project
        self.branch = ''
        self.src_path = ''
        if hasattr(app, 'config'):
            self.project = app.config.srclink_project
            self.branch = app.config.srclink_branch
            self.src_path = app.config.srclink_src_path

        (self.project_path,
         self.project_url) = self.get_project_vars_from_urlstr(self.project)

        self.app = app
        self.pagename = pagename
        self.templatename = templatename
        self.context = context
        self.doctree = doctree
        self.path = path

    @staticmethod
    def get_context_class_for_url(urlstr):
        clsmap = Context.build_type_context_clsmap()
        _url = urlparse.urlparse(urlstr)
        for _, Context_cls in clsmap.items():
            if Context_cls.url_matches(url=_url):
                return Context_cls

    def get_project_vars_from_urlstr(self, urlstr):
        _url = urlparse.urlparse(urlstr)
        _path = _url.path
        if _path.startswith('git@') or _path.startswith('hg@'):
            _path = _url.path.lstrip("git@").lstrip("hg@")
        pathcomp = _path.split('/')
        if not _url.netloc:
            if len(pathcomp) == 2:
                project_path = '/'.join(pathcomp)
            elif len(pathcomp) == 3:
                project_path = '/'.join(pathcomp[-2:])
            else:
                raise Exception("Unable to parse %r (%r)" % (urlstr, _path))
        else:
            project_path = '/'.join([p for p in _path.split('/') if p][:2])  # user/project
        if self.type_:
            project_url = 'https://{type_}/{project_path}'.format(
                type_=self.type_,
                project_path=project_path)
        else:
            project_url = project_path
        return (project_path, project_url)

    @staticmethod
    def get_srclink_context(urlstr):
        if '/' in urlstr or '@' in urlstr:
            klass = Context.get_context_class_for_url(urlstr)
            if klass is None:
                if '/' in urlstr:
                    klass = GitHubContext
        else:
            klass = GitHubContext
        if klass is None:
            raise Exception("No Context could be found for %r" % urlstr)
        return klass

    lookup = get_srclink_context

    @staticmethod
    def url_matches(urlstr=None, url=None):
        return

    def get_url(self, view=None, path=None, conf=None, tmpl=None, rev=None):
        """
        This is the general pattern of both bitbucket and github,
        with ``s/(blame|annotate)/annotate/`` and ``s/(blob|src)/src/``
        e.g. for testing.

        subclasses can/should generate an ``https://url`` ::

            class Projectsite(Context):
                def get_url(app, view, path):
                    return tmpl.format(**conf)

        """
        if tmpl is None:
            tmpl = self.src_view_url_tmpl
        if rev is None:
            rev = self.branch
        if conf is None:
            # Sphinx-gallery specific
            if "auto" in path and "plot" in path:
                path = path.replace(".rst", ".py")
                path = path.replace("_examples_auto", "examples")
                path = path.replace("_tutorials_auto", "tutorials")

            if "installation" in path:
                path = path.replace("doc/installation.rst", "INSTALLATION.rst")

            conf = dict(
                project_url=self.project_url,
                project=self.project,
                project_path=self.project_path,
                view=view,
                rev=rev,
                src_path=self.src_path,
                path=path,)

        return tmpl.format(**conf)

    @classmethod
    def get_context(
            cls,
            app=None,
            pagename=None,
            templatename=None,
            context=None,
            doctree=None,
            path=None):
        if app and doctree:
            path = os.path.relpath(doctree.get('source'), app.builder.srcdir)
        ctxt = cls(app, pagename, templatename, context, doctree, path)
        context.update(ctxt.to_dict())
        return context

    def get_srclinks(self, project=None):
        if project is None:
            project = self.project
        conf = {
            'project': self.project,
            'project_path': self.project_path,
            'project_url': self.project_url,
            'src_type': self.src_type,
            'rev': self.branch,
        }
        c = collections.OrderedDict()
        c['srclink_project_path'] = self.project_path

        if self.src_url_tmpl is not None:
            c['srclink_src_url'] = self.src_url_tmpl.format(**conf)

        c['srclink_src_rev'] = conf['rev']
        if self.src_rev_url_tmpl is not None:
            c['srclink_src_rev_url'] = self.src_rev_url_tmpl.format(**conf)

        if self.src_https_tmpl is not None:
            c['srclink_src_https_url'] = self.src_https_tmpl.format(**conf)
        if self.src_ssh_tmpl is not None:
            c['srclink_src_ssh_url'] = self.src_ssh_tmpl.format(**conf)
        if self.src_native_tmpl is not None:
            c['srclink_src_native_url'] = self.src_native_tmpl.format(**conf)
        if self.src_clonecmd_tmpl is not None:
            c['srclink_src_clone_cmd'] = self.src_clonecmd_tmpl.format(
                **conf)
        return c

    def get_page_action_urls(self, app=None, path=None):
        """
        Subclases can/should override this method
        """
        c = collections.OrderedDict()
        c['show_srclink_url'] = self.get_url('blob', path)
        c['edit_srclink_url'] = self.get_url('edit', path)
        c['history_srclink_url'] = self.get_url('commits', path)
        c['annotate_srclink_url'] = self.get_url('blame', path)
        return c

    def to_dict(self):
        c = collections.OrderedDict()
        c.update(self.get_srclinks(self.project))
        c.update(self.get_page_action_urls(app=self.app, path=self.path))
        return c

    @staticmethod
    def build_type_context_clsmap(default_context=None):
        clsmap = collections.OrderedDict([
            ((cls.type_, cls.src_type), cls) for cls in [
                BitbucketGitSourcelinkContext,
                BitbucketSourcelinkContext,
                GitHubContext,
                ]])
        if default_context is None:
            default_context = Context
        clsmap[None] = default_context
        return clsmap


class GitHubContext(Context):
    type_ = 'github.com'
    src_type = 'git'
    src_https_tmpl = 'https://github.com/{project_path}'
    src_url_tmpl = src_https_tmpl
    src_rev_url_tmpl = 'https://github.com/{project_path}/tree/{rev}'
    src_ssh_tmpl = 'ssh://git@github.com/{project_path}'
    src_native_tmpl = 'git://github.com/{project_path}'
    src_view_url_tmpl = (
        'https://github.com/{project_path}/{view}/{rev}/{src_path}{path}')

    def match_issue_regex(textstr):
        raise NotImplementedError()
        issue_regexes = [
            r"#[\d]+",  # TODO: {,10}
            r"GH[\d]+",
        ]
        for rgx in issue_regexes:
            yield re.match(rgx, textstr)

    @staticmethod
    def is_github_url(urlstr=None, url=None):
        if url is None:
            url = urlparse.urlparse(urlstr)
        if (url.netloc in ('git@github.com', 'github.com')
                or url.path.startswith('git@github.com/')):
            return True

    url_matches = is_github_url


class BitbucketSourcelinkContext(Context):
    type_ = 'bitbucket.org'

    src_type = 'hg'
    src_https_tmpl = 'https://bitbucket.org/{project_path}'
    src_url_tmpl = src_https_tmpl
    src_rev_url_tmpl = 'https://bitbucket.org/{project_path}/src/{rev}'
    #src_ssh_git_tmpl = 'ssh://git@bitbucket.org/{project_path}'
    #src_src_git_tmpl = 'git://bitbucket.org/{project_path}'
    src_ssh_hg_tmpl = 'ssh://hg@bitbucket.org/{project_path}'
    src_ssh_tmpl = src_ssh_hg_tmpl
    src_src_hg_tmpl = 'hg://bitbucket.org/{project_path}'
    src_native_tmpl = src_src_hg_tmpl

    src_view_url_tmpl = (
        'https://bitbucket.org/{project_path}/{view}/{rev}/{src_path}{path}')

    def __init__(self, *args, **kwargs):
        super(BitbucketSourcelinkContext, self).__init__(*args, **kwargs)
        _url = urlparse.urlparse(self.project_url)
        if (_url.path.startswith('hg@bitbucket.org/')
                or _url.netloc == 'hg@bitbucket.org'):
            self.src_type = 'hg'
            self.src_ssh_tmpl = self.src_ssh_hg_tmpl
            self.src_native_tmpl = self.src_src_hg_tmpl
        elif (_url.path.startswith('git@bitbucket.org/')
              or _url.netloc == 'git@bitbucket.org'):
            self.src_type = 'git'
            self.src_ssh_tmpl = self.src_ssh_git_tmpl
            self.src_native_tmpl = self.src_src_git_tmpl

    @staticmethod
    def is_bitbucket_url(urlstr=None, url=None):
        if url is None:
            url = urlparse.urlparse(urlstr)
        if 'bitbucket.org' in url.netloc:
            return True

    url_matches = is_bitbucket_url

    def get_page_action_urls(self, app=None, path=None):
        c = collections.OrderedDict()
        c['show_srclink_url'] = self.get_url('blob', path)
        c['edit_srclink_url'] = self.get_url('edit', path)
        c['history_srclink_url'] = self.get_url('commits', path)
        c['annotate_srclink_url'] = self.get_url('blame', path)
        return c


class BitbucketGitSourcelinkContext(BitbucketSourcelinkContext):
    type_ = 'bitbucket.org'
    src_type = 'git'
    src_https_tmpl = 'https://bitbucket.org/{project_path}'
    src_url_tmpl = src_https_tmpl
    src_ssh_tmpl = 'ssh://git@bitbucket.org/{project_path}'
    src_native_tmpl = 'git://bitbucket.org/{project_path}'

    src_view_url_tmpl = (
        'https://bitbucket.org/{project_path}/{view}/{rev}/{src_path}{path}')

    @staticmethod
    def is_bitbucket_git_url(urlstr=None, url=None):
        if url is None:
            url = urlparse.urlparse(urlstr)
        if ('git@bitbucket.org' == url.netloc
                or url.path.startswith('git@bitbucket.org/')):
            return True

    url_matches = is_bitbucket_git_url


def html_page_context(app=None,
                      pagename=None,
                      templatename=None,
                      context=None,
                      doctree=None,
                      urlstr=None):
    """
    Sphinx HTML page context
    """
    if templatename != 'page.html':
        return
    if not doctree:
        # warnings.warn("doctree is None")
        return
    if app:
        if urlstr is None:
            urlstr = app.config.srclink_project
    if not urlstr:
        warnings.warn("srclink_project not specified")
        return
    ctxt = Context.lookup(urlstr)
    variables = ctxt.get_context(app, pagename, templatename, context, doctree)
    # import pprint
    # pprint.pprint(variables)
    context.update(variables)


def setup(app):
    app.add_config_value('srclink_project', '', True)
    app.add_config_value('srclink_branch', 'master', True)
    app.add_config_value('srclink_src_path', '', True)  # 'eg' "docs/"
    app.connect('html-page-context', html_page_context)

import unittest


class TestRepositorycontexts(unittest.TestCase):

    def test_srclinks(self):
        class App(object):

            def __init__(self, **kwargs):
                self.__dict__.update(kwargs)

        def build_srclink_context(project=None, branch=None, src_path=None):
            if project is None:
                project = "PROJECT"
            if branch is None:
                branch = "BRANCH"
            if src_path is None:
                src_path = "docs/"
            context = {
                'srclink_project': project,
                'srclink_branch': branch,
                'srclink_src_path': src_path,
            }
            return context

        def test_html_page_context(project,
                                   project_branch,
                                   project_src_path,
                                   project_path,
                                   project_url=None,
                                   project_https_url=None,
                                   project_ssh_url=None,
                                   project_native_url=None):
            if project_url is None:
                project_url = project
            context = build_srclink_context(
                project=project,
                branch=project_branch,
                src_path=project_src_path)
            app = App(config=App(**context), builder=App(srcdir='srcdir'))
            pagename = 'TEST'
            templatename = 'page.html'
            doctree = {'source': 'DOCTREE_SOURCE'}
            html_page_context(
                app,
                pagename,
                templatename,
                context,
                doctree)
            self.assertTrue(context)
            self.assertEqual(context['srclink_project'], project)
            self.assertEqual(context['srclink_branch'], project_branch)
            self.assertEqual(context['srclink_src_path'], project_src_path)
            self.assertEqual(context['srclink_project_path'], project_path)
            self.assertEqual(context['srclink_src_url'], project_url)
            self.assertEqual(
                context['srclink_src_https_url'],
                project_https_url)
            self.assertEqual(context['srclink_src_ssh_url'], project_ssh_url)
            self.assertEqual(
                context['srclink_src_native_url'],
                project_native_url)
            return context

        test_html_page_context(
            "github.com/westurner/dotfiles",
            "develop",
            "doc/",
            "westurner/dotfiles",
            project_url="https://github.com/westurner/dotfiles",
            project_https_url="https://github.com/westurner/dotfiles",
            project_ssh_url="ssh://git@github.com/westurner/dotfiles",
            project_native_url="git://github.com/westurner/dotfiles")

        test_html_page_context(
            "westurner/dotfiles",
            "develop",
            "doc/",
            "westurner/dotfiles",
            project_url="https://github.com/westurner/dotfiles",
            project_https_url="https://github.com/westurner/dotfiles",
            project_ssh_url="ssh://git@github.com/westurner/dotfiles",
            project_native_url="git://github.com/westurner/dotfiles")

        test_html_page_context(
            "https://github.com/westurner/dotfiles",
            "develop",
            "doc/",
            "westurner/dotfiles",
            project_https_url="https://github.com/westurner/dotfiles",
            project_ssh_url="ssh://git@github.com/westurner/dotfiles",
            project_native_url="git://github.com/westurner/dotfiles")

        test_html_page_context(
            "https://bitbucket.org/westurner/dotfiles",
            "default",
            "doc/",
            "westurner/dotfiles",
            project_https_url="https://bitbucket.org/westurner/dotfiles",
            project_ssh_url="ssh://hg@bitbucket.org/westurner/dotfiles",
            project_native_url="hg://bitbucket.org/westurner/dotfiles")

        test_html_page_context(
            "git@bitbucket.org/westurner/dotfiles",
            "master",
            "doc/",
            "westurner/dotfiles",
            project_url="https://bitbucket.org/westurner/dotfiles",
            project_https_url="https://bitbucket.org/westurner/dotfiles",
            project_ssh_url="ssh://git@bitbucket.org/westurner/dotfiles",
            project_native_url="git://bitbucket.org/westurner/dotfiles")


def main(argv=None):
    args = argv
    if args is None:
        args = sys.argv

    loglevel = logging.INFO
    if '-v' in args:
        args.remove('-v')
        loglevel = logging.DEBUG
    elif '-q' in args:
        args.remove('-q')
        loglevel = logging.ERROR
    logging.basicConfig(level=loglevel)

    log = logging.getLogger(__name__)
    log.info(args)

    if '-t' in args:
        args.remove('-t')
        logging.basicConfig(level=logging.DEBUG)
        sys.exit(unittest.main(argv=args))

    fields = ['srclink_project', 'srclink_branch', 'srclink_src_path']
    fieldtuples = itertools.izip_longest(fields, args[1:])
    conf = collections.OrderedDict(fieldtuples)
    print("### CONF ###")
    print(json.dumps(conf, indent=4))

    urlstr = conf.get('srclink_project')

    print("### context ###")
    output = html_page_context(urlstr=urlstr, context=conf)
    print(output)

    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
