import re
from numpy import int32, int64

def parse_index(string, shift=0):
    if not string:
        return []
    assert re.match(r'(\d+(\-\d+)?)(,\d+(\-\d+)?)*', string)
    segs = string.strip().split(',')
    indices = []
    for seg in segs:
        m = re.match(r'(\d+)\-(\d+)', seg)
        if m:
            indices += list(range(int(m.group(1)), int(m.group(2)) + 1))
        else:
            indices.append(int(seg))
    return [i+shift for i in indices]

def rev_parse_index(l):
    templ = list(l)
    for ele in templ:
        assert type(ele) in (int, int32, int64)
    templ.sort()
    start = templ.pop(0)
    end = start
    s = ''
    while templ:
        ele = templ.pop(0)
        # print ele, start, end
        if ele == end:
            continue
        if ele == end + 1:
            end = ele
        else:
            if start == end:
                s = s + '%d,' % start
            else:
                s = s + '%d-%d,' % (start, end)
            start = ele
            end = start
    if start == end:
        s = s + '%d' % start
    else:
        s = s + '%d-%d' % (start, end)
    return s

def test():
    while True:
        s = raw_input()
        if not s:
            break
        print parse_index(s)

def test2():
    # while True:
        # s = raw_input()
        # if not s:
            # break
        # print rev_parse_index(parse_index(s))
    
    import random
    for it in range(100000):
        l = [random.randint(1,100) for _ in range(10)]
        # print l
        try:
            assert parse_index(rev_parse_index(l)) == sorted(list(set(l)))
        except AssertionError:
            print l
            print rev_parse_index(l)
            print parse_index(rev_parse_index(l))
            print sorted(list(set(l)))
            break

if __name__ == '__main__':
    from codetester import tester
    # test()
    tester(test2)

