import filecmp


def filediff(f1, f2):
    """
    Returns a string summarizing differences between two text files.
    Returns the empty string if no differences.

    Args:
        f1: path to first file
        f2: path to second file
    """
    msg = []
    if not filecmp.cmp(f1, f2):
        msg = [f'{f1} differs from expected:']
        found = False
        with open(f1) as f1, open(f2) as f2:
            for n, (line1, line2) in enumerate(zip(f1, f2), 1):
                if line1 != line2:
                    found = True
                    msg.append(f'first difference at line {n}:')
                    msg.append(f' - {line2}')
                    msg.append(f' + {line1}')
            if not found:
                f1.seek(0)
                f2.seek(0)
                ll1, ll2 = len(f1.readlines()), len(f2.readlines())
                assert ll1 != ll2
                msg.append(
                    f'first {n} lines match but {f1} has '
                    f'{ll1} lines and {f2} has {ll2}')

    return '\n'.join(msg)

