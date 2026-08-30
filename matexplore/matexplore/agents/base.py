# -*- coding: utf-8 -*-
"""智能体基类：提供 run(ctx) 协议 + 可插拔 LLM 接口。

LLM provider=none 时走确定性推理(各子类实现)；配置了 provider 时，
可调用 _ask_llm() 让大模型对数据结果做自然语言推理(假设润色/候选点评)。
"""
import os


class Agent:
    name = "base"

    def __init__(self, cfg):
        self.cfg = cfg
        self.llm = None
        if cfg.get("llm", {}).get("provider", "none") != "none":
            self.llm = cfg.llm

    def run(self, ctx):
        raise NotImplementedError

    def _ask_llm(self, prompt, system=None):
        """调用外部 LLM(OpenAI 兼容)。未配置时抛异常，由调用方回退。"""
        if self.llm is None:
            raise RuntimeError("LLM not configured")
        import urllib.request, json
        key = os.environ.get(self.llm.api_key_env, "")
        body = {"model": self.llm.model, "messages": []}
        if system:
            body["messages"].append({"role": "system", "content": system})
        body["messages"].append({"role": "user", "content": prompt})
        req = urllib.request.Request(self.llm.base_url + "/chat/completions",
                                     data=json.dumps(body).encode(),
                                     headers={"Content-Type": "application/json",
                                              "Authorization": "Bearer " + key})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = json.loads(r.read())
        return data["choices"][0]["message"]["content"]

