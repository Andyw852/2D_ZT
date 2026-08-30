# -*- coding: utf-8 -*-
"""配置加载。

为避免对 PyYAML 的硬依赖(目标计算环境可能没有)，内置一个仅覆盖本框架所需
YAML 子集的解析器；若环境有 PyYAML 则优先使用。
"""
import os


class DotDict(dict):
    def __getattr__(self, k):
        try:
            return self[k]
        except KeyError:
            raise AttributeError(k)


def _to_dot(v):
    """载入时一次性深转：dict->DotDict, list->list(DotDict)，之后所有访问返回同一对象，mutation 可持久化。"""
    if isinstance(v, dict):
        return DotDict({k: _to_dot(x) for k, x in v.items()})
    if isinstance(v, list):
        return [_to_dot(x) for x in v]
    return v


def _strip_comment(line):
    in_s = None
    for i, ch in enumerate(line):
        if ch in ('"', "'"):
            if in_s is None:
                in_s = ch
            elif in_s == ch:
                in_s = None
        elif ch == '#' and in_s is None:
            return line[:i]
    return line


def _parse_scalar(s):
    s = s.strip()
    if s == '':
        return None
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s[1:-1]
    if len(s) >= 2 and s[0] == "'" and s[-1] == "'":
        return s[1:-1]
    low = s.lower()
    if low == 'true':
        return True
    if low == 'false':
        return False
    if low in ('null', '~', 'none'):
        return None
    if s.startswith('[') and s.endswith(']'):
        inner = s[1:-1].strip()
        if not inner:
            return []
        return [_parse_scalar(x) for x in _split_top(inner, ',')]
    try:
        if ('.' in s) or ('e' in low) or ('E' in s):
            return float(s)
        return int(s)
    except ValueError:
        return s


def _split_top(s, sep):
    out, depth, cur = [], 0, ''
    for ch in s:
        if ch in '[{(':
            depth += 1
        elif ch in ']})':
            depth -= 1
        if ch == sep and depth == 0:
            out.append(cur)
            cur = ''
        else:
            cur += ch
    out.append(cur)
    return out


def _parse(text):
    raw = []
    for ln in text.splitlines():
        ln = _strip_comment(ln)
        if ln.strip() == '':
            continue
        indent = len(ln) - len(ln.lstrip(' '))
        raw.append((indent, ln.strip()))

    def parse_block(i, indent):
        result = {}
        while i < len(raw):
            ind, content = raw[i]
            if ind < indent:
                break
            if ind > indent:
                i += 1
                continue
            if content.startswith('- '):
                return parse_list(i, indent)
            if ':' in content:
                key, val = content.split(':', 1)
                key = key.strip().strip('"').strip("'")
                val = val.strip()
                if val == '':
                    if i + 1 < len(raw) and raw[i + 1][0] > indent:
                        nind = raw[i + 1][0]
                        if raw[i + 1][1].startswith('- '):
                            sub, i = parse_list(i + 1, nind)
                        else:
                            sub, i = parse_block(i + 1, nind)
                        result[key] = sub
                    else:
                        result[key] = None
                        i += 1
                else:
                    result[key] = _parse_scalar(val)
                    i += 1
            else:
                i += 1
        return result, i

    def parse_list(i, indent):
        items = []
        while i < len(raw):
            ind, content = raw[i]
            if ind < indent or ind > indent:
                break
            if content.startswith('- '):
                item = content[2:].strip()
                if item == '':
                    sub, i = parse_block(i + 1, indent + 2)
                    items.append(sub)
                else:
                    items.append(_parse_scalar(item))
                    i += 1
            else:
                break
        return items, i

    result, _ = parse_block(0, 0)
    return result


def load_config(path=None):
    if path is None:
        path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "config.yaml")
    with open(path, "r") as f:
        text = f.read()
    try:
        import yaml
        data = yaml.safe_load(text)
    except ImportError:
        data = _parse(text)
    return _to_dot(data)

