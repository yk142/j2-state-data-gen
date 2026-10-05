function X = j2Sobol(n, d, seed)
%J2SOBOL Sobol 低食い違い列（最大 4 次元）を自前実装し、[0,1)^d の n 点を返す。
%   X = J2SOBOL(n, d)         標準の Sobol 列（先頭点は 0）
%   X = J2SOBOL(n, d, seed)   seed を与えると、次元ごとのランダム XOR シフト（digital shift）を
%                             かけて再現可能にばらつかせる（低食い違い性は保たれる）。
%
%   Statistics and Machine Learning Toolbox の sobolset は使用禁止のため自前で実装する。
%   方向数は Joe–Kuo の初期値（次元 1: van der Corput、次元 2: s=1, 次元 3: s=2, a=1,
%   m=[1 3]、次元 4: s=3, a=1, m=[1 3 1]）。点生成は Gray code 法。
if nargin < 3, seed = []; end
assert(d >= 1 && d <= 4, 'j2Sobol:dim', '対応次元は 1〜4 です');
B = 32;
S  = [0 1 2 3];   A = [0 0 1 1];
M0 = {[], 1, [1 3], [1 3 1]};
V = zeros(d, B, 'uint32');
for k = 1:B, V(1,k) = bitshift(uint32(1), B-k); end
for j = 2:d
    s = S(j);  a = A(j);  m = M0{j};
    for k = 1:s, V(j,k) = bitshift(uint32(m(k)), B-k); end
    for k = s+1:B
        v = bitxor(V(j,k-s), bitshift(V(j,k-s), -s));
        for l = 1:s-1
            if bitget(a, s-l), v = bitxor(v, V(j,k-l)); end
        end
        V(j,k) = v;
    end
end
shift = zeros(1, d, 'uint32');
if ~isempty(seed)
    rs = RandStream('twister','Seed',seed);
    shift = uint32(floor(rand(rs,1,d) * 2^32));
end
x = zeros(1, d, 'uint32');
X = zeros(n, d);
X(1,:) = double(bitxor(x, shift)) / 2^32;
for i = 2:n
    c = find(~bitget(uint32(i-2), 1:B), 1);        % (i-2) の最下位の 0 ビット位置
    x = bitxor(x, V(:,c)');
    X(i,:) = double(bitxor(x, shift)) / 2^32;
end
end
