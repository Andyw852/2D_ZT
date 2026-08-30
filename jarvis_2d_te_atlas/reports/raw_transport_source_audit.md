# Phase G1：JARVIS 原始输运数据源审计

## jarvis-tools 官方下载机制

- 源码：jarvis/db/figshare.py
- get_db_info()["dft_2d"] = 以下 4 项：
  1. URL：https://ndownloader.figshare.com/files/38521268
  2. ZIP 内部 JSON：d2-12-12-2022.json
  3. 描述：Obtaining 2D dataset 1.1k ...
  4. Reference：https://doi.org/10.1016/j.commatsci.2025.114063
- 本地缓存文件名：get_request_data() 里 zfile = js_tag + ".zip"，即 d2-12-12-2022.json.zip，
  默认存放于 get_cache_dir("jarvis_data") = $HOME/.cache/atomgptlab/jarvis_data（可被 ATOMGPTLAB_CACHE 覆盖）。

## FigShare 状态

- 本执行环境对 figshare.com / api.figshare.com / ndownloader.figshare.com 全部返回 HTTP 403
  （AWS ALB 网络层阻断），已在 Phase 1 reports/download_probe.json 记录。
- 因此 data("dft_2d") 无法在本环境下载 figshare JSON。

## 本机缓存搜索（Phase G2）

- $HOME/.cache/atomgptlab/jarvis_data：空（不存在 d2-12-12-2022.json.zip）
- $HOME/.jarvis：不存在
- site-packages / /tmp / 项目目录：未发现 d2-12-12-2022.json 或对应 ZIP
- 结论：本机无缓存，需要从官方接口重新获取。

## 关键发现：官方 NIST 静态 XML 服务器包含原始 3 本征值（Phase G3）

- OPTIMADE 每条记录 _jarvis_reference 字段指向：
  https://www.ctcms.nist.gov/~knc6/static/JARVIS-DFT/{JID}
- 该地址返回完整 JARVIS-DFT XML 记录（HTTP 200，官方 NIST 静态文件，Cloudflare 可达）。
- XML 中 <boltztrap_info> 字段包含 p/n 的 seeb / cond / pf / kappa 各 3 个本征值（逗号分隔字符串）：
    pseeb  = 251.58,251.58,119.37
    pcond  = 19412.37,19412.32,0.51
    ppf    = 1228.62,1228.62,0.01
    pkappa = 91235212000000.0,91234998000000.0,542820440.0
    nseeb  = -147.13,-147.13,-11.33
    ncond  = 44252.96,44252.86,0.32
    npf    = 957.91,957.91,0.0
    nkappa = 109512420000000.0,109512169999999.98,436661230.0
- XML 中 <effective_mass> 字段包含 electron/hole 有效质量 3 个本征值：
    electron_mass_300K = 0.51,0.51,7358.37
    hole_mass_300K = 1.78,1.78,20893.98

## OPTIMADE scalar 是否隐藏数组字段（Phase G4）

- 检查 OPTIMADE structures 记录全部 attributes：仅结构类字段为 list（lattice_vectors、
  cartesian_site_positions、species 等），_jarvis_efg / structure_features 为空 list。
- 输运字段（n/p-Seebeck、ncond/pcond、nkappa/pkappa、n-powerfact/p-powerfact）全部为 scalar float。
- 结论：OPTIMADE 不提供原始数组，数组在 NIST 静态 XML 服务器中。

## 2D 材料 3 本征值的物理结构（以 MoS2 为例）

- pseeb = [251.58, 251.58, 119.37]：两个面内（in-plane）本征值几乎简并，一个面外（out-of-plane-like）本征值不同。
- pcond = [19412.37, 19412.32, 0.51]：面内两个很大且相等，面外约 0（符合二维材料跨真空层不导电）。
- 注意：数据库只保存本征值、不保存本征向量，因此不能确定哪个值对应 x/y/z 晶向，
  只能识别「两个简并面内 + 一个面外」这种 2D 特征谱型（且不能默认最小值为 z，见任务第 39 节）。
