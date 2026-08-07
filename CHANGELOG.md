# Changelog

## [0.9.0](https://github.com/fg-labs/bwa-mem3/compare/v0.8.0...v0.9.0) (2026-08-06)


### Features

* **mem:** make --compat=bwa-mem selectable ([#360](https://github.com/fg-labs/bwa-mem3/issues/360)) ([22d604d](https://github.com/fg-labs/bwa-mem3/commit/22d604d98a9d55d982dbb5e121bb96c9ef3f61b6))


### Bug Fixes

* **bam:** propagate bam_aux_append failures instead of dropping tags ([#367](https://github.com/fg-labs/bwa-mem3/issues/367)) ([9e7c600](https://github.com/fg-labs/bwa-mem3/commit/9e7c6001b866ec3f044c8f54a5f63d0c6af28c60))
* **bam:** render pa:f from the same definition as the SAM writer ([#366](https://github.com/fg-labs/bwa-mem3/issues/366)) ([fd7e3fe](https://github.com/fg-labs/bwa-mem3/commit/fd7e3fe04a5345d68bd60d73095602f3493c8848)), closes [#365](https://github.com/fg-labs/bwa-mem3/issues/365)
* **pair:** derive FLAG 0x2 from a[0] by default, matching both upstreams ([#363](https://github.com/fg-labs/bwa-mem3/issues/363)) ([a9a96cc](https://github.com/fg-labs/bwa-mem3/commit/a9a96ccfaf4cd7e6f74aa6a510c4ecf97ab3013e))


### Documentation

* **changelog:** backtick SAM header tags to stop @-mentioning users ([#359](https://github.com/fg-labs/bwa-mem3/issues/359)) ([769baec](https://github.com/fg-labs/bwa-mem3/commit/769baeccd5341eb2b30f3f83be5c22b4186bf023))
* **readme:** state bwa-mem2 parity and surface --compat ([#364](https://github.com/fg-labs/bwa-mem3/issues/364)) ([668b21e](https://github.com/fg-labs/bwa-mem3/commit/668b21ecfa2645f9a37eba63dd6a9f7e093794dd))

## [0.8.0](https://github.com/fg-labs/bwa-mem3/compare/v0.7.0...v0.8.0) (2026-08-03)


### ⚠ BREAKING CHANGES

* **meth:** `bwa-mem3 mem --meth` now writes SAM text instead of BAM. Add `--bam` to any script that depended on the old behavior. Pipelines that pipe into samtools need no change — samtools autodetects SAM text.

### Features

* **bam:** warn that in-process compressed BAM is single-threaded ([#287](https://github.com/fg-labs/bwa-mem3/issues/287)) ([9cd414d](https://github.com/fg-labs/bwa-mem3/commit/9cd414dcfd4f98b75f1fb71f570a769143455923))
* **mem:** --compat target enum with record + header byte-identity to bwa-mem2 ([#277](https://github.com/fg-labs/bwa-mem3/issues/277)) ([265b427](https://github.com/fg-labs/bwa-mem3/commit/265b42769731b09d8cca1221148c5f2be98faeb6))
* **mem:** band mate rescue to a k-mer anchor diagonal (--rescue-kmer) ([#335](https://github.com/fg-labs/bwa-mem3/issues/335)) ([6ce2b70](https://github.com/fg-labs/bwa-mem3/commit/6ce2b70509446c874fb7f6bd037bd6553e0a1e23))
* **meth:** add --meth-tags to select which Bismark tags are emitted ([#333](https://github.com/fg-labs/bwa-mem3/issues/333)) ([3ec4500](https://github.com/fg-labs/bwa-mem3/commit/3ec4500fc3e26c9e3140db1c0d4aa81fd84d1f8d)), closes [#331](https://github.com/fg-labs/bwa-mem3/issues/331)
* **meth:** let --bam choose the output container under --meth ([#341](https://github.com/fg-labs/bwa-mem3/issues/341)) ([ed4400c](https://github.com/fg-labs/bwa-mem3/commit/ed4400c5f09302160366403789eddbcf7d41a5f5))


### Bug Fixes

* **bntseq:** don't poison the pac-fetch buffer in release builds ([#263](https://github.com/fg-labs/bwa-mem3/issues/263)) ([1252f1b](https://github.com/fg-labs/bwa-mem3/commit/1252f1bf4245fa4f2245c64157109e6cc17bd88e))
* **bsw:** compute the 8-bit band clamp in wide arithmetic ([#270](https://github.com/fg-labs/bwa-mem3/issues/270)) ([00ade6d](https://github.com/fg-labs/bwa-mem3/commit/00ade6d9930bf4d5bf0609dd3480b58b4fee7509))
* **bsw:** correct 8-bit banded-SW z-drop and seed clamp at high seed scores (all tiers) ([#273](https://github.com/fg-labs/bwa-mem3/issues/273)) ([e722ed0](https://github.com/fg-labs/bwa-mem3/commit/e722ed06dd71d0186c401ce2f0e163d32b0a7ad7))
* **build:** generate header dependencies instead of a hand-maintained list ([#299](https://github.com/fg-labs/bwa-mem3/issues/299)) ([414ac7e](https://github.com/fg-labs/bwa-mem3/commit/414ac7e614db17a3564ae5f9cfdab4f5e6d54ad4))
* **chain:** keep sizeof(mem_chain_t) at 48 B to preserve bwa-mem2 chaining parity ([#268](https://github.com/fg-labs/bwa-mem3/issues/268)) ([a9d9445](https://github.com/fg-labs/bwa-mem3/commit/a9d9445cb984115887184e5bdb0a5ddc41e74708))
* check debug_macro_flag_lint's directory argument before cd ([#352](https://github.com/fg-labs/bwa-mem3/issues/352)) ([14d8096](https://github.com/fg-labs/bwa-mem3/commit/14d8096fec6e657e6008044267354f57c42c3fea))
* **coderabbit:** list main in auto_review.base_branches ([#358](https://github.com/fg-labs/bwa-mem3/issues/358)) ([63add93](https://github.com/fg-labs/bwa-mem3/commit/63add93a2a2cc2f02f216ff4ddbf1cb07ede634d))
* **header:** emit one default `@HD` on every output path ([#291](https://github.com/fg-labs/bwa-mem3/issues/291)) ([56889e1](https://github.com/fg-labs/bwa-mem3/commit/56889e1f5733b72f333fa12faf0ea82157f92d3f)), closes [#288](https://github.com/fg-labs/bwa-mem3/issues/288)
* **index:** cap each pread() at 1GiB so large indexes load on macOS ([#259](https://github.com/fg-labs/bwa-mem3/issues/259)) ([1db6697](https://github.com/fg-labs/bwa-mem3/commit/1db6697cb5456a5905c1e7ce736f4f4d63431334))
* **index:** resolve the auto memory budget from the host, not a 32 GiB cap ([#300](https://github.com/fg-labs/bwa-mem3/issues/300)) ([a8ac368](https://github.com/fg-labs/bwa-mem3/commit/a8ac3683c521bd7f388bb201f75b9a23af1280d4))
* **kswv:** zero query padding in the NEON and AVX-512BW 8-bit mate-rescue kernels ([#290](https://github.com/fg-labs/bwa-mem3/issues/290)) ([b8d5aa3](https://github.com/fg-labs/bwa-mem3/commit/b8d5aa38856d7149bcdf9161a1d5ca38ec5a5bdc))
* **mem:** guard allocations with xassert so OOM checks survive NDEBUG ([#312](https://github.com/fg-labs/bwa-mem3/issues/312)) ([8898782](https://github.com/fg-labs/bwa-mem3/commit/88987823179e8aa179a53d15b8743fa17f22b947))
* **mem:** make the batch-size cap opt-in so default batching matches bwa-mem2 ([#298](https://github.com/fg-labs/bwa-mem3/issues/298)) ([5a85305](https://github.com/fg-labs/bwa-mem3/commit/5a85305fee874ba16aa3e47d10eedbf4e7ea248e))
* **meth:** derive NM/MD from the scoring matrix so conversions are not counted ([#332](https://github.com/fg-labs/bwa-mem3/issues/332)) ([6aa4675](https://github.com/fg-labs/bwa-mem3/commit/6aa4675636673d94fdbbe768671e4ff56e4cad07)), closes [#327](https://github.com/fg-labs/bwa-mem3/issues/327)
* **meth:** diagnose orphaned --meth/--meth-tags values instead of aligning to them ([#334](https://github.com/fg-labs/bwa-mem3/issues/334)) ([a18b19b](https://github.com/fg-labs/bwa-mem3/commit/a18b19bf48a45d1bfb303645b04e73f458365743)), closes [#331](https://github.com/fg-labs/bwa-mem3/issues/331)
* **meth:** emit MQ:i and HN:i from the --meth BAM writer ([#304](https://github.com/fg-labs/bwa-mem3/issues/304)) ([70d1270](https://github.com/fg-labs/bwa-mem3/commit/70d1270666b3d50d5a51c94b795a7570cd28531a)), closes [#296](https://github.com/fg-labs/bwa-mem3/issues/296)
* **meth:** reject mate rescues whose SW alignment runs past the read ([#258](https://github.com/fg-labs/bwa-mem3/issues/258)) ([6ecc61b](https://github.com/fg-labs/bwa-mem3/commit/6ecc61b30a3a38c225cd7641de4ae7ed372e25fa))
* **sam:** make the SAM-A9 rid check opt-in so the recompute is skipped ([#330](https://github.com/fg-labs/bwa-mem3/issues/330)) ([87f3bd6](https://github.com/fg-labs/bwa-mem3/commit/87f3bd6b56eac5b01b97e496fd5730106f52a6b0))
* **simd:** map _MM_HINT_T0 to L1, not L3, on arm64 ([#262](https://github.com/fg-labs/bwa-mem3/issues/262)) ([34cb389](https://github.com/fg-labs/bwa-mem3/commit/34cb3891d5afd773279900b8428d0b5bf956a921))


### Performance

* **bntseq:** O(1) contig bucket table for bns_pos2rid ([#275](https://github.com/fg-labs/bwa-mem3/issues/275)) ([e0291ec](https://github.com/fg-labs/bwa-mem3/commit/e0291ec3a352e79b3dcd16a254f4ea6ce451313c))
* **bntseq:** unpack the 2-bit reference with a byte-&gt;4-base LUT ([#274](https://github.com/fg-labs/bwa-mem3/issues/274)) ([a9b570f](https://github.com/fg-labs/bwa-mem3/commit/a9b570fb302b5b365612164b85ddad6339b1b28c))
* **bsw:** admit high-h0 pairs to the 8-bit banded-SW tier ([#321](https://github.com/fg-labs/bwa-mem3/issues/321)) ([35f3238](https://github.com/fg-labs/bwa-mem3/commit/35f3238c1f4d52d2b0e94fc037db68a22d6314f6))
* **bsw:** drop a redundant compare in the 8-bit row argmax ([#266](https://github.com/fg-labs/bwa-mem3/issues/266)) ([ed34d5d](https://github.com/fg-labs/bwa-mem3/commit/ed34d5d9787313aab2ed012dcce1a8ec04870f76))
* **bsw:** fuse the SBT pre-pass into the 8-bit DP loop with two LUTs ([#280](https://github.com/fg-labs/bwa-mem3/issues/280)) ([94cbdf2](https://github.com/fg-labs/bwa-mem3/commit/94cbdf271a08c70ebf20abe645b412df33bc8ff5))
* **bsw:** pack the int8 lane groups by max(len1,len2) ([#283](https://github.com/fg-labs/bwa-mem3/issues/283)) ([7898254](https://github.com/fg-labs/bwa-mem3/commit/78982548c6b306daaea945a9333a70390ea34fb2))
* **chain:** hoist chain bounds and memoize the per-chain log() ([#265](https://github.com/fg-labs/bwa-mem3/issues/265)) ([7703634](https://github.com/fg-labs/bwa-mem3/commit/7703634884c5cc7c70f1479b121c2782ce527932))
* **chain:** skip side-effect-free comparisons in mem_chain_flt ([#326](https://github.com/fg-labs/bwa-mem3/issues/326)) ([6ca3c7e](https://github.com/fg-labs/bwa-mem3/commit/6ca3c7ef6bea30f09f927867bcf145c8cd20b9de))
* **io:** pool per-read string fields in a per-chunk arena ([#293](https://github.com/fg-labs/bwa-mem3/issues/293)) ([c63eba0](https://github.com/fg-labs/bwa-mem3/commit/c63eba080791c864e95000aa979b9733cfb7ec02))
* **kswv:** recover the query end after the row, via block checkpoints ([#328](https://github.com/fg-labs/bwa-mem3/issues/328)) ([f099e79](https://github.com/fg-labs/bwa-mem3/commit/f099e79d1c4e0cf2f06b00fe23ee60f013f20ba2))
* **kswv:** stop testing for query padding on columns that cannot have any ([#324](https://github.com/fg-labs/bwa-mem3/issues/324)) ([90be8df](https://github.com/fg-labs/bwa-mem3/commit/90be8df08bb2823b95dfa27897da92c08c17e7d9))
* **main:** calibrate proc_freq without sleeping a full second ([#295](https://github.com/fg-labs/bwa-mem3/issues/295)) ([052a84a](https://github.com/fg-labs/bwa-mem3/commit/052a84a20f229a915592f975181de9eaf4bf7ff2))
* **mem:** fuse the seed and extend kt_for passes ([#264](https://github.com/fg-labs/bwa-mem3/issues/264)) ([d4127e9](https://github.com/fg-labs/bwa-mem3/commit/d4127e9cb7a02b8b033402d622954d183c71d656))
* **mem:** read a pestat cohort in slices so compute starts sooner ([#305](https://github.com/fg-labs/bwa-mem3/issues/305)) ([212ed34](https://github.com/fg-labs/bwa-mem3/commit/212ed34d3a0254db78e7a7db3c029582f53e65d7))
* **mem:** reserve the cohort accumulator from task_size instead of doubling ([#342](https://github.com/fg-labs/bwa-mem3/issues/342)) ([80e40ad](https://github.com/fg-labs/bwa-mem3/commit/80e40ad351c41a1a23e4c7a8278cdee809e57b56))
* **mem:** size chaining scratch per thread instead of per read ([#297](https://github.com/fg-labs/bwa-mem3/issues/297)) ([3d66be4](https://github.com/fg-labs/bwa-mem3/commit/3d66be4da53acea114cf45de25c347c798e70f6a))
* **mem:** size the cohort ramp so it cannot outrun the reader ([#315](https://github.com/fg-labs/bwa-mem3/issues/315)) ([160adae](https://github.com/fg-labs/bwa-mem3/commit/160adae9b418adcf320247a016bb9caea3a186f0))
* **mem:** size the per-chunk read pool from the actual read count ([#279](https://github.com/fg-labs/bwa-mem3/issues/279)) ([3848339](https://github.com/fg-labs/bwa-mem3/commit/3848339b55706bcece1fc18ccbc884ed27651448))
* **meth:** share the generic writer's per-record scratch with the meth writer ([#306](https://github.com/fg-labs/bwa-mem3/issues/306)) ([8fade15](https://github.com/fg-labs/bwa-mem3/commit/8fade150e18195fd863d7bdda9b07a07684b5607))
* **rescue:** drop redundant mate-query copies in batched mate-SW ([#284](https://github.com/fg-labs/bwa-mem3/issues/284)) ([32f17e5](https://github.com/fg-labs/bwa-mem3/commit/32f17e5809fa744a0d74d20d24365d45b893d8bb))
* **rescue:** single-position anchor index, plus an opt-in --rescue-skip gate ([#349](https://github.com/fg-labs/bwa-mem3/issues/349)) ([e7e00af](https://github.com/fg-labs/bwa-mem3/commit/e7e00af4471a6f400a34346cf868016223f5e8a8))
* **sam,seed:** micro-cleanups in SAM formatting and seeding ([#286](https://github.com/fg-labs/bwa-mem3/issues/286)) ([82257ab](https://github.com/fg-labs/bwa-mem3/commit/82257abcf33a6056c49490e444165d48ee803973))
* **seed:** correct the FM-index checkpoint prefetch targets ([#267](https://github.com/fg-labs/bwa-mem3/issues/267)) ([6455cad](https://github.com/fg-labs/bwa-mem3/commit/6455cadea07e60859740b3174ebb79fbfd397ed2))
* **seed:** reuse sortSMEMs counting-sort scratch across batches ([#278](https://github.com/fg-labs/bwa-mem3/issues/278)) ([394f8f8](https://github.com/fg-labs/bwa-mem3/commit/394f8f8110f7d15be7ef2ca38c335590aa1e0284))
* **sort:** move the pdqsort/total-order dedup sort behind --fast ([#257](https://github.com/fg-labs/bwa-mem3/issues/257)) ([34b3446](https://github.com/fg-labs/bwa-mem3/commit/34b34466bce9412842d7eed1da3c772d38bad87d))
* **sort:** run pdqsort in the dedup sorts when the comparator sees no tie ([#261](https://github.com/fg-labs/bwa-mem3/issues/261)) ([7481c4f](https://github.com/fg-labs/bwa-mem3/commit/7481c4ffa7ac9a1a55b249fd0f7269b5e4c38f9b))


### Refactoring

* **bsw:** derive extension gaps from M again (bwa-mem2 compatibility) ([#256](https://github.com/fg-labs/bwa-mem3/issues/256)) ([771376e](https://github.com/fg-labs/bwa-mem3/commit/771376e3839275910704972d9a1a1915766b1eec))
* **bsw:** drop the unreachable 8-bit re-baselining in the AVX2/AVX-512 kernels ([#272](https://github.com/fg-labs/bwa-mem3/issues/272)) ([f1a0da6](https://github.com/fg-labs/bwa-mem3/commit/f1a0da65bab8f1cdca2192b9e56671767df0cabd))
* **bsw:** drop the unreachable 8-bit score re-baselining ([#271](https://github.com/fg-labs/bwa-mem3/issues/271)) ([ecebf2b](https://github.com/fg-labs/bwa-mem3/commit/ecebf2b34fb26bdddd212007bf58b25cf893032e))
* **chain:** fuse mem_chain_weight's two seed sweeps into one ([#322](https://github.com/fg-labs/bwa-mem3/issues/322)) ([500c9cb](https://github.com/fg-labs/bwa-mem3/commit/500c9cbcb3a62c023673c58fdc5cc985619dc4f8))
* **chain:** remove dead work from the chaining path ([#285](https://github.com/fg-labs/bwa-mem3/issues/285)) ([8de7e5b](https://github.com/fg-labs/bwa-mem3/commit/8de7e5bb91229b7b6efb9684d82f0cb6d22f4c4f))
* **chain:** tidy mem_chain_flt's filter loops ([#323](https://github.com/fg-labs/bwa-mem3/issues/323)) ([095b401](https://github.com/fg-labs/bwa-mem3/commit/095b401be04aa65348b5f65da7438f0347e3edc1))
* **header:** one generated-`@SQ` builder, one scanner, one record iterator ([#292](https://github.com/fg-labs/bwa-mem3/issues/292)) ([37e5e4a](https://github.com/fg-labs/bwa-mem3/commit/37e5e4a6a18c66ed3e4be5e818a4b57f5dafff86))
* **meth:** name the native-regen decision and document the hypothesis invariant ([#337](https://github.com/fg-labs/bwa-mem3/issues/337)) ([97eb3be](https://github.com/fg-labs/bwa-mem3/commit/97eb3be3a92c3f5963b41e1ad0cfef6a5aeaa55c))


### Documentation

* correct the developer guide's account of the source-only lints ([#354](https://github.com/fg-labs/bwa-mem3/issues/354)) ([a8825f5](https://github.com/fg-labs/bwa-mem3/commit/a8825f52a59a8abe4eed9673c1c3e77a33c046ff))
* **equivalence:** reflect restored bwa-mem2 byte-identity on the drop-in profile ([#276](https://github.com/fg-labs/bwa-mem3/issues/276)) ([8dbf218](https://github.com/fg-labs/bwa-mem3/commit/8dbf218adeeb1215e04619d921bbec18e0746d6f))
* **kswv:** the phase-1 boundary-mask cost is measured now, not unmeasurable ([#308](https://github.com/fg-labs/bwa-mem3/issues/308)) ([136cc9f](https://github.com/fg-labs/bwa-mem3/commit/136cc9ffac6ce857933b9bd7850549402e63de47))
* **memory:** measure the -t batch multiplier cost and correct the resident index figure ([#319](https://github.com/fg-labs/bwa-mem3/issues/319)) ([804f20b](https://github.com/fg-labs/bwa-mem3/commit/804f20bd41b806af4873325583681a0bc1bad1dc))
* **test:** describe the source-only lints' input contract ([#339](https://github.com/fg-labs/bwa-mem3/issues/339)) ([de77830](https://github.com/fg-labs/bwa-mem3/commit/de77830153e019b829d8f33ccd609e55703e149e))
* **test:** replace the dead phiX example for running a regression test ([#345](https://github.com/fg-labs/bwa-mem3/issues/345)) ([d3d9a53](https://github.com/fg-labs/bwa-mem3/commit/d3d9a5322123fc45320e017ae52650af19a94cf2))

## [0.7.0](https://github.com/fg-labs/bwa-mem3/compare/v0.6.0...v0.7.0) (2026-07-23)


### Features

* **bsw:** add BWAMEM3_DUMP_PAIRS to observe extension pairs ([#239](https://github.com/fg-labs/bwa-mem3/issues/239)) ([c4d7480](https://github.com/fg-labs/bwa-mem3/commit/c4d748084e4d844e1448f3c7f7b0c5039f9d8406))
* **meth:** add --meth=taps for TAPS methylation chemistry ([#228](https://github.com/fg-labs/bwa-mem3/issues/228)) ([0b8163d](https://github.com/fg-labs/bwa-mem3/commit/0b8163d2632bd168d2848617858176ac9ce5c127))
* **meth:** drop cross-chemistry seeds by read number ([#225](https://github.com/fg-labs/bwa-mem3/issues/225)) ([6e21a37](https://github.com/fg-labs/bwa-mem3/commit/6e21a37e40e0013916863673622b01fd694dbf58))
* **meth:** NEUTRAL --meth-scoring (TAPS default) + batched kswv rescue ([#243](https://github.com/fg-labs/bwa-mem3/issues/243)) ([f36dc5a](https://github.com/fg-labs/bwa-mem3/commit/f36dc5a05e820985ed3f9f802f5a88840c3ea978))


### Bug Fixes

* **ext:** tighten ungapped fast-path mismatch bound to the strict form ([#231](https://github.com/fg-labs/bwa-mem3/issues/231)) ([5c53540](https://github.com/fg-labs/bwa-mem3/commit/5c5354042b2036882fa8840b3ad62a78f7eed853))
* **fast:** enable --extend-mate-concordant and raise chain cap to 20 for non-meth --fast ([#218](https://github.com/fg-labs/bwa-mem3/issues/218)) ([50a3293](https://github.com/fg-labs/bwa-mem3/commit/50a329390753e77640449eafa564ecf13ade6ec3))
* **mem:** guard extension-buffer growth against int32 offset overflow ([#237](https://github.com/fg-labs/bwa-mem3/issues/237)) ([6acd6b4](https://github.com/fg-labs/bwa-mem3/commit/6acd6b43d4d1712b71f9244f85addc3ca8cd920b))
* **mem:** size sa_coord realloc in elements, not SMEM count ([#235](https://github.com/fg-labs/bwa-mem3/issues/235)) ([7fa501c](https://github.com/fg-labs/bwa-mem3/commit/7fa501c1a6e67d419b1a216f5c8b282385d518c6))
* **meth:** recover repeat placement via both-hypothesis mate rescue + PE MAPQ hardening ([#219](https://github.com/fg-labs/bwa-mem3/issues/219)) ([6705bbf](https://github.com/fg-labs/bwa-mem3/commit/6705bbfa89c2ab40f4fb81f2f6a6a9518c1cf35c))
* **meth:** scale --meth scoring defaults with the match score ([#227](https://github.com/fg-labs/bwa-mem3/issues/227)) ([9811767](https://github.com/fg-labs/bwa-mem3/commit/98117672a6bd0780a9a39ef8b4e50b4f1a6a539e))


### Performance

* **bandedSWA:** fewer ops in the AVX2/AVX-512 extend cores (x86-only) ([#255](https://github.com/fg-labs/bwa-mem3/issues/255)) ([d52842f](https://github.com/fg-labs/bwa-mem3/commit/d52842feef5b8c7cb6137e552c02c2234ddd206c))
* **bsw:** 128-bit (SSE/NEON) 16-bit prepass via byte-LUT + sign-extend ([#223](https://github.com/fg-labs/bwa-mem3/issues/223)) ([239a334](https://github.com/fg-labs/bwa-mem3/commit/239a334b0d6dcbc60ee4e76076353ea79e342479))
* **bsw:** AVX-512 16-bit prepass via 32-entry permutexvar LUT ([#222](https://github.com/fg-labs/bwa-mem3/issues/222)) ([0e8c0b6](https://github.com/fg-labs/bwa-mem3/commit/0e8c0b6674229243c26196e421caeb48580f2c23))
* **index:** parallel pread of the FM-index for faster startup ([#248](https://github.com/fg-labs/bwa-mem3/issues/248)) ([8d88fe5](https://github.com/fg-labs/bwa-mem3/commit/8d88fe5700be4bbbb757976b51a03bcd7e8b2032))
* **kswv:** fold the freed-cell pair into a per-row active_frread on the u16 kernels ([#254](https://github.com/fg-labs/bwa-mem3/issues/254)) ([eaf649f](https://github.com/fg-labs/bwa-mem3/commit/eaf649fcc5943d91f9d70bd4e215a6d489871bc5))
* **kswv:** hoist the freed-cell target into a per-row active_frread vector ([#252](https://github.com/fg-labs/bwa-mem3/issues/252)) ([fc0eb91](https://github.com/fg-labs/bwa-mem3/commit/fc0eb9132179ccac9d79ffe6c480951160582e60))
* **kswv:** hoist the u8 boundary test per-row on NEON and AVX-512 non-meth ([#253](https://github.com/fg-labs/bwa-mem3/issues/253)) ([a25e3b7](https://github.com/fg-labs/bwa-mem3/commit/a25e3b7bf7b3444cc0eadccce6243cd23279e076))
* **kswv:** skip the ref-boundary vbsl on all-real rows in 8-bit mate rescue (+7.5%) ([#244](https://github.com/fg-labs/bwa-mem3/issues/244)) ([78420a5](https://github.com/fg-labs/bwa-mem3/commit/78420a5553b85a7a3fbb98f61017354c8bf7e20a))
* **mem:** resolve SA lookups across reads, not per read ([#236](https://github.com/fg-labs/bwa-mem3/issues/236)) ([79a2425](https://github.com/fg-labs/bwa-mem3/commit/79a24254168c4cbef086dcfeb073ecbd5632e438))
* **meth:** score mate rescue with one read#-derived matrix, not both ([#224](https://github.com/fg-labs/bwa-mem3/issues/224)) ([1baa39b](https://github.com/fg-labs/bwa-mem3/commit/1baa39b17c69584263694fda5ddb666bcc65ad71))
* **meth:** vectorize both-hypothesis mate rescue on the batched kswv path ([#221](https://github.com/fg-labs/bwa-mem3/issues/221)) ([46aeaae](https://github.com/fg-labs/bwa-mem3/commit/46aeaaeb818ae8e5272f29ba8807e76db59f7d5b)), closes [#219](https://github.com/fg-labs/bwa-mem3/issues/219)
* **seed:** prefetch the correct cp_occ lines in the bwtSeed third round ([#242](https://github.com/fg-labs/bwa-mem3/issues/242)) ([c5e3827](https://github.com/fg-labs/bwa-mem3/commit/c5e3827df83b5dd85a998423d03cb93bec99263f))


### Documentation

* **best-practices:** add anti-patterns for marginally-mappable input ([#232](https://github.com/fg-labs/bwa-mem3/issues/232)) ([1344ab8](https://github.com/fg-labs/bwa-mem3/commit/1344ab86826cf6b3ba4a563bb4b482addb506428))
* **ext:** record the measured evidence behind the ungapped fast-path bound ([#233](https://github.com/fg-labs/bwa-mem3/issues/233)) ([878bda0](https://github.com/fg-labs/bwa-mem3/commit/878bda0972509f44aa2073bafd0fb51997002b33))
* **mem:** scope long-read claims to what actually works ([#240](https://github.com/fg-labs/bwa-mem3/issues/240)) ([ac090bc](https://github.com/fg-labs/bwa-mem3/commit/ac090bc13a3b2e37217b0989e204fce15f9fd07e))
* **meth:** correct --meth override claims and document -T 40 split filtering ([#226](https://github.com/fg-labs/bwa-mem3/issues/226)) ([50a8d3f](https://github.com/fg-labs/bwa-mem3/commit/50a8d3f7fb56a631ec030eef4edfd1042619710f))
* **profiles:** recommend --skip-contained-ext for short-read pipelines ([#241](https://github.com/fg-labs/bwa-mem3/issues/241)) ([37fa409](https://github.com/fg-labs/bwa-mem3/commit/37fa4098e5180a8a0b99e878b65095562ec3ec15))

## [0.6.0](https://github.com/fg-labs/bwa-mem3/compare/v0.5.0...v0.6.0) (2026-07-16)


### Features

* **version:** report whether mimalloc is the active allocator ([#217](https://github.com/fg-labs/bwa-mem3/issues/217)) ([2cd9fe9](https://github.com/fg-labs/bwa-mem3/commit/2cd9fe9c27185af98d2f6a0b59c04431a50eea0a))


### Performance

* **bam:** cut per-record work in the --bam writer (mem_aln_to_bam) ([#212](https://github.com/fg-labs/bwa-mem3/issues/212)) ([837f40e](https://github.com/fg-labs/bwa-mem3/commit/837f40e1e974217ee1c45f7a4ceea6b3b6c00821))
* **bsw:** drop dead per-row H1/H2 setup stores in the SW batch wrappers ([#211](https://github.com/fg-labs/bwa-mem3/issues/211)) ([039b09f](https://github.com/fg-labs/bwa-mem3/commit/039b09f8ecd11e724b098a1729d42703c2ce058d))
* **dedup:** drop provably-dead exact-duplicate passes in mem_sort_dedup_patch ([#205](https://github.com/fg-labs/bwa-mem3/issues/205)) ([c544258](https://github.com/fg-labs/bwa-mem3/commit/c5442582e44da347d8a81e9b58a70a18ed0467a4))
* **fmi:** arm64 lockstep for third-pass reseeding (bwtSeedStrategy) ([#215](https://github.com/fg-labs/bwa-mem3/issues/215)) ([c5f69ec](https://github.com/fg-labs/bwa-mem3/commit/c5f69ec0aaa8d2fae290228c7b24ee774f740d62))
* **kswv:** remove write-only Hmax scratch buffer from batched SW kernels ([#214](https://github.com/fg-labs/bwa-mem3/issues/214)) ([b14ef49](https://github.com/fg-labs/bwa-mem3/commit/b14ef4993ee9a60813b0412b8a5d5d40aa19b66b))
* **mem:** add unique-mapper fast paths to mem_chain_flt and mem_mark_primary_se ([#209](https://github.com/fg-labs/bwa-mem3/issues/209)) ([1700aab](https://github.com/fg-labs/bwa-mem3/commit/1700aabc263a8776dc12b66acbb566bda67a37ba))
* **mem:** drop redundant per-call scratch allocations in mem_gen_alt ([#216](https://github.com/fg-labs/bwa-mem3/issues/216)) ([3418525](https://github.com/fg-labs/bwa-mem3/commit/3418525750b15df76dd0d5c5d3634dab4c8f61d4))
* **mem:** pool per-read scratch in get_sa_entries_prefetch and mem_reg2aln ([#208](https://github.com/fg-labs/bwa-mem3/issues/208)) ([5747a13](https://github.com/fg-labs/bwa-mem3/commit/5747a131cde9c22490c5ace470b0903d395ad130))
* **mem:** remove dead/instrumentation work from the extension hot path ([#206](https://github.com/fg-labs/bwa-mem3/issues/206)) ([439be97](https://github.com/fg-labs/bwa-mem3/commit/439be9709844a0557ad01ec7cdf3b542e9acf66d))
* **seed:** replace sortSMEMs qsort with a counting sort by rid ([#207](https://github.com/fg-labs/bwa-mem3/issues/207)) ([4356b5c](https://github.com/fg-labs/bwa-mem3/commit/4356b5ce50348a2421b2b8e2517d12f642cfb8e8))
* **seed:** vectorize backwardExt occ-counting on arm64 (NEON) + hoist invariants ([#210](https://github.com/fg-labs/bwa-mem3/issues/210)) ([64a55a2](https://github.com/fg-labs/bwa-mem3/commit/64a55a25c18350a66a93438d64f3d93a77259a8a))


### Refactoring

* **mem:** static cleanup — remove dead code, silence int64 format warnings ([#213](https://github.com/fg-labs/bwa-mem3/issues/213)) ([8c57b7e](https://github.com/fg-labs/bwa-mem3/commit/8c57b7ef6b0ff45d3b6d1b17031c9f0c46879199))


### Documentation

* retitle LICENSE for BWA-MEM3 and add Fulcrum copyright ([#203](https://github.com/fg-labs/bwa-mem3/issues/203)) ([089a956](https://github.com/fg-labs/bwa-mem3/commit/089a956375b27372d812a757cc87663526126fc8))

## [0.5.0](https://github.com/fg-labs/bwa-mem3/compare/v0.4.0...v0.5.0) (2026-07-04)


### Features

* add opt-in --seed-order seed reordering (default off, byte-identical) ([#186](https://github.com/fg-labs/bwa-mem3/issues/186)) ([04749a1](https://github.com/fg-labs/bwa-mem3/commit/04749a16fb2e2a4f0bd8a318fc6bdd6439e2689f))
* add opt-in --smem-dedup (dedup identical SMEMs before chaining) ([#187](https://github.com/fg-labs/bwa-mem3/issues/187)) ([1384972](https://github.com/fg-labs/bwa-mem3/commit/1384972417677943b68b8292181017223ad27876))
* **mem:** add --adaptive-band (chain-geometry adaptive banding) for long reads ([#194](https://github.com/fg-labs/bwa-mem3/issues/194)) ([4fe92a6](https://github.com/fg-labs/bwa-mem3/commit/4fe92a6ea55ebe8617cabed79902bceeabecbd2c))
* **mem:** add --extend-mate-concordant; fix --fast --meth placement regression ([#195](https://github.com/fg-labs/bwa-mem3/issues/195)) ([c9ffef1](https://github.com/fg-labs/bwa-mem3/commit/c9ffef1551d13ef22b104a75c5be785857e01efa))
* **mem:** add --fast speed preset ([#189](https://github.com/fg-labs/bwa-mem3/issues/189)) ([a946af8](https://github.com/fg-labs/bwa-mem3/commit/a946af8fcd0309faed952662123f88d152b3b585))
* **mem:** add --max-extend-chains and bundle it into --fast ([#193](https://github.com/fg-labs/bwa-mem3/issues/193)) ([e39b3d4](https://github.com/fg-labs/bwa-mem3/commit/e39b3d497e3e10535b2d9eabb9bb7fa48d1c1379))
* **mem:** add --skip-contained-ext and enable it under --fast ([#192](https://github.com/fg-labs/bwa-mem3/issues/192)) ([2d2b2b4](https://github.com/fg-labs/bwa-mem3/commit/2d2b2b437b70cdc5c445ea83da7cb158af654cd1))


### Bug Fixes

* **bandedSWA:** 8-bit SW drops query-end gscore/gtle on zero-score-row exit ([#198](https://github.com/fg-labs/bwa-mem3/issues/198)) ([611e21b](https://github.com/fg-labs/bwa-mem3/commit/611e21b924d09ccbd7e2e3d18092a90dc798e4b2))
* **bandedSWA:** getScores{8,16} must not scribble padding past numPairs ([#199](https://github.com/fg-labs/bwa-mem3/issues/199)) ([9aae808](https://github.com/fg-labs/bwa-mem3/commit/9aae808c0a721bf88dfc4c486c8bf8b785cdd508))


### Performance

* **bandedSWA:** gate the getScores overshoot guard to sub-slice callers ([#201](https://github.com/fg-labs/bwa-mem3/issues/201)) ([162e909](https://github.com/fg-labs/bwa-mem3/commit/162e90923488fea4634c00a990156794940b9797))


### Refactoring

* **kswv:** drop duplicate F warm-up prefetch in kswv512_16 ([#191](https://github.com/fg-labs/bwa-mem3/issues/191)) ([6e4cf2b](https://github.com/fg-labs/bwa-mem3/commit/6e4cf2b34a527309e8dc675e941a9fd548feadf1))


### Documentation

* **changelog:** backfill the 0.4.0 breaking-change notice ([#183](https://github.com/fg-labs/bwa-mem3/issues/183)) ([e1c381a](https://github.com/fg-labs/bwa-mem3/commit/e1c381ad4b27e7e3bcf26dbbfeea130cff6e2af0))
* **changelog:** render the live changelog, not the frozen NEWS.md ([#184](https://github.com/fg-labs/bwa-mem3/issues/184)) ([12a46d8](https://github.com/fg-labs/bwa-mem3/commit/12a46d8317fd3dabb031f835348ebdf126e76601))
* **contributing:** document breaking-change commit footers ([#181](https://github.com/fg-labs/bwa-mem3/issues/181)) ([4ebd122](https://github.com/fg-labs/bwa-mem3/commit/4ebd1220ffda0a0da544bfe0a9bab94840d73a8c))
* **release:** describe the release-please flow, not manual tagging ([#182](https://github.com/fg-labs/bwa-mem3/issues/182)) ([fc0c4b6](https://github.com/fg-labs/bwa-mem3/commit/fc0c4b65f356e6571c9f33f563a13d9fc6bd3f3e))

## [0.4.0](https://github.com/fg-labs/bwa-mem3/compare/v0.3.0...v0.4.0) (2026-06-27)


### ⚠ BREAKING CHANGES

* **index:** `bwa-mem3 index` no longer writes the unpacked `.0123` reference file by default ([#177](https://github.com/fg-labs/bwa-mem3/pull/177)). `bwa-mem3 mem` now reconstructs reference bases from the packed `.pac` on demand ("pac-fetch") and ignores any `.0123` present, so the file is redundant for bwa-mem3 itself. External tools that read `.0123` directly — most notably sharing a single index with **bwa-mem2** — will break, since the expected file is now absent. To restore the old on-disk layout, re-run indexing with the new opt-in flag `index --emit-unpacked-ref`. Alignment output is byte-for-byte identical; the change is purely to the index artifact set.


### Features

* **mem:** --min-ext-len opt-in filter to skip extension of short seeds ([#169](https://github.com/fg-labs/bwa-mem3/issues/169)) ([13db252](https://github.com/fg-labs/bwa-mem3/commit/13db252ad4c7c39eff18bd3674430aac51570346))
* **meth:** native bisulfite (BS-seq) alignment via --meth (D3) ([#174](https://github.com/fg-labs/bwa-mem3/issues/174)) ([a0296b1](https://github.com/fg-labs/bwa-mem3/commit/a0296b1766f17aab85b720fef1cf7f3984735105))


### Performance

* **mem:** pac-fetch the reference from .pac instead of loading/building .0123 ([#177](https://github.com/fg-labs/bwa-mem3/issues/177)) ([9c4bbf2](https://github.com/fg-labs/bwa-mem3/commit/9c4bbf217d3cbf02143a2fbd9fb75329126d9010))
* **meth:** batched (SIMD) asymmetric mate rescue (closes [#173](https://github.com/fg-labs/bwa-mem3/issues/173)) ([#175](https://github.com/fg-labs/bwa-mem3/issues/175)) ([f146a18](https://github.com/fg-labs/bwa-mem3/commit/f146a18a32232a6ea05b040b10917a464a2d264c))


### Documentation

* **book:** document the libdeflate build prerequisite (incl. AL2023) ([#172](https://github.com/fg-labs/bwa-mem3/issues/172)) ([c2b6ec7](https://github.com/fg-labs/bwa-mem3/commit/c2b6ec78bc8fadb1aa63f99aee4f90671f840e7d))
* **book:** recommend -y 0 (drop 3rd-round seeding) as an opt-in speed knob ([#171](https://github.com/fg-labs/bwa-mem3/issues/171)) ([ca9ac1f](https://github.com/fg-labs/bwa-mem3/commit/ca9ac1f2289c200708193c3f6ebbe094f230d79d))
* collapse the mdbook sidebar into nested, foldable sections ([#180](https://github.com/fg-labs/bwa-mem3/issues/180)) ([c6ac47b](https://github.com/fg-labs/bwa-mem3/commit/c6ac47b9ebc3146144377f4a56e689b3fdf11946))
* deep mdbook cleanup — dedup, consolidate, and tighten ([#179](https://github.com/fg-labs/bwa-mem3/issues/179)) ([54d6d11](https://github.com/fg-labs/bwa-mem3/commit/54d6d1161c374d7eb7fcef0eaac7c1156fe0bfb4))
* **meth:** disclose collapsed-mode placement drift vs bwameth.py ([#178](https://github.com/fg-labs/bwa-mem3/issues/178)) ([96f29e2](https://github.com/fg-labs/bwa-mem3/commit/96f29e257e9acd9072caa246b1589c8943babc26))
* **settings-profiles:** note repeat-aggregating downstream caveat for -m 10 ([#168](https://github.com/fg-labs/bwa-mem3/issues/168)) ([38fd1ec](https://github.com/fg-labs/bwa-mem3/commit/38fd1ec9717e29cace9eac33bd643a7da001da4e))

## [0.3.0](https://github.com/fg-labs/bwa-mem3/compare/v0.2.2...v0.3.0) (2026-06-21)


### Features

* **bsw:** make the 8-bit h0-prefix seed unsigned [0,255] ([#151](https://github.com/fg-labs/bwa-mem3/issues/151)) ([9f51c5f](https://github.com/fg-labs/bwa-mem3/commit/9f51c5f521612a7c3e1c0dae9d71db38be0121b8))
* **bsw:** recover the 8-bit banded Smith–Waterman path for reads ≥128 bp ([#140](https://github.com/fg-labs/bwa-mem3/issues/140)) ([155a916](https://github.com/fg-labs/bwa-mem3/commit/155a91632de602663377197ff1168d200a24f344))
* **kswv:** AVX2 16-bit mate-rescue kernel (kswv256_16) ([#162](https://github.com/fg-labs/bwa-mem3/issues/162)) ([9107b82](https://github.com/fg-labs/bwa-mem3/commit/9107b825189da87822bccac152e4a98daf00aac1))
* **meth:** carry original-reference `@SQ` `M5`/`UR` and `@CO`/`@PG` into `--meth` headers ([#139](https://github.com/fg-labs/bwa-mem3/issues/139)) ([e94ad8b](https://github.com/fg-labs/bwa-mem3/commit/e94ad8b847b6637b95aea6dff3e07f560e22a532))
* **prof:** off-by-default --profile stage-timing instrumentation ([#152](https://github.com/fg-labs/bwa-mem3/issues/152)) ([83cf7ab](https://github.com/fg-labs/bwa-mem3/commit/83cf7ab9d8e4589354e16654726fbdfaa1b4b63c))
* **reader:** content-detecting FASTQ reader fast path (libdeflate BGZF) ([#128](https://github.com/fg-labs/bwa-mem3/issues/128)) ([cdd71bf](https://github.com/fg-labs/bwa-mem3/commit/cdd71bf4590cf2e1a5d908d6e852f780c93c3a8a))


### Bug Fixes

* **bsw:** bound getScores8/16 prefetch reads to the padding contract ([#150](https://github.com/fg-labs/bwa-mem3/issues/150)) ([87ed5d4](https://github.com/fg-labs/bwa-mem3/commit/87ed5d4a78b6b78789bca45be5649373576ee668))
* **fmi:** widen mem_lim to int64 and guard SA-entry allocations ([#156](https://github.com/fg-labs/bwa-mem3/issues/156)) ([2d18c1e](https://github.com/fg-labs/bwa-mem3/commit/2d18c1ed00de5e4a60fa31d15fec04c7b74877a0))
* **kthread:** drive kt_for with a persistent worker pool ([#154](https://github.com/fg-labs/bwa-mem3/issues/154)) ([26b24e7](https://github.com/fg-labs/bwa-mem3/commit/26b24e7bac9fd7f73921b902ba98a628776087b6))
* **meth:** emit `-R` read group as `@RG` header in `--meth` mode ([#137](https://github.com/fg-labs/bwa-mem3/issues/137)) ([ccd1fc5](https://github.com/fg-labs/bwa-mem3/commit/ccd1fc56d26bfd54ab1f9604845332353afd75a6))
* **seeding:** widen SMEM read positions from int16_t to int32_t ([#142](https://github.com/fg-labs/bwa-mem3/issues/142)) ([037c418](https://github.com/fg-labs/bwa-mem3/commit/037c418b3b3c49deeb41bd57e37eacc7895e6fca))
* **test:** make meth layer-2 FAIL diagnostics reachable under set -e ([#133](https://github.com/fg-labs/bwa-mem3/issues/133)) ([d2d6688](https://github.com/fg-labs/bwa-mem3/commit/d2d66880bce736900b81e7610cbae1d1556ba7ea))


### Performance

* **bsw:** AVX2 SIMD tuning for the Smith-Waterman kernels ([#161](https://github.com/fg-labs/bwa-mem3/issues/161)) ([458b216](https://github.com/fg-labs/bwa-mem3/commit/458b21653cb6f84808032360a85a78f291064a4e))
* **bsw:** NEON SIMD tuning for the Smith-Waterman kernels ([#160](https://github.com/fg-labs/bwa-mem3/issues/160)) ([d971ff0](https://github.com/fg-labs/bwa-mem3/commit/d971ff0694d8a90cf836dfff873534c9d3a160a4))
* **bsw:** prefetch next batch's ref/query in the AVX2 8-bit wrapper ([#163](https://github.com/fg-labs/bwa-mem3/issues/163)) ([e6082a0](https://github.com/fg-labs/bwa-mem3/commit/e6082a05fc6af07da26dae90210ea651b601dd91))
* **bsw:** short-circuit the inert per-row re-baseline scan ([#147](https://github.com/fg-labs/bwa-mem3/issues/147)) ([8e284e0](https://github.com/fg-labs/bwa-mem3/commit/8e284e0cf274e8c77c78457e372c1ab92de69f3f))
* **bsw:** vectorize the per-row epilogue side-channel loop ([#149](https://github.com/fg-labs/bwa-mem3/issues/149)) ([403aeb7](https://github.com/fg-labs/bwa-mem3/commit/403aeb76ad8ec2664736998c018a7e033af0c070))
* **fmi:** size SA-entry staging buffers to the exact write count ([#157](https://github.com/fg-labs/bwa-mem3/issues/157)) ([aa0fe33](https://github.com/fg-labs/bwa-mem3/commit/aa0fe3335c58badf97893fae32f67ff02f69217a))
* **read:** vendored zlib-ng inflate + chunk cap + 3rd pipeline worker ([#153](https://github.com/fg-labs/bwa-mem3/issues/153)) ([5cf89e3](https://github.com/fg-labs/bwa-mem3/commit/5cf89e3c60e01bec089183e511efa27327414035))
* **sw:** reassociate affine-gap recurrences on NEON (kswv + bandedSWA) ([#166](https://github.com/fg-labs/bwa-mem3/issues/166)) ([a02fcb4](https://github.com/fg-labs/bwa-mem3/commit/a02fcb446574d5b5d03abdbf73c9b129deead2d4))


### Refactoring

* **bsw:** derive extension gaps from H (standard Gotoh), not M ([#141](https://github.com/fg-labs/bwa-mem3/issues/141)) ([f715fbd](https://github.com/fg-labs/bwa-mem3/commit/f715fbd92b4c4acabe7b7ab8dde35d5ad7fe211c))
* **bsw:** drop the dead qlen[] parameter from the 8-bit kernels ([#143](https://github.com/fg-labs/bwa-mem3/issues/143)) ([42321df](https://github.com/fg-labs/bwa-mem3/commit/42321df6100f04451b8a4d5ed7f214ea28522ec8))
* **bsw:** remove dead SW code paths (SORT_PAIRS, non-CORE macros, SSE2 polyfill) ([#148](https://github.com/fg-labs/bwa-mem3/issues/148)) ([ad8937e](https://github.com/fg-labs/bwa-mem3/commit/ad8937ef05fb48992ca61be70c21d4b6e245a7b3))


### Documentation

* add memory budgeting and data-type tuning guide ([#145](https://github.com/fg-labs/bwa-mem3/issues/145)) ([7127d80](https://github.com/fg-labs/bwa-mem3/commit/7127d80f7204c7951bf74d8191a883c4e0f00a82))
* add situational --supp-rep-hard-cap 20 note for SV-aware pipelines ([#134](https://github.com/fg-labs/bwa-mem3/issues/134)) ([e20bc0e](https://github.com/fg-labs/bwa-mem3/commit/e20bc0e28d3501a732e9b1eba88a9920ec877bae))
* **perf:** refresh reference-architecture table to v0.3.0 (a02fcb4) ([#167](https://github.com/fg-labs/bwa-mem3/issues/167)) ([2cd7bfa](https://github.com/fg-labs/bwa-mem3/commit/2cd7bfa1cfb5392856905df3b88cebeda2a1d659))
* **perf:** what drives the speedup, full perf-PR catalog, and fix the stale RTD build ([#155](https://github.com/fg-labs/bwa-mem3/issues/155)) ([0b01fb7](https://github.com/fg-labs/bwa-mem3/commit/0b01fb73a2e9a8b9e8899f30b2afd94c9e6a1418))
* recommend -s 0 for --meth Pass-2 re-seeding in settings profiles ([#132](https://github.com/fg-labs/bwa-mem3/issues/132)) ([45e02e0](https://github.com/fg-labs/bwa-mem3/commit/45e02e0095f862c99eaf99273b2912ff2229ec9b))
* recommend a recent compiler on aarch64, with measured NEON numbers ([#165](https://github.com/fg-labs/bwa-mem3/issues/165)) ([b591684](https://github.com/fg-labs/bwa-mem3/commit/b59168479ddb0fddba86cd0eb0e61ee0a1c8fcb0))
* settings profiles (bwa drop-in vs recommended) ([#131](https://github.com/fg-labs/bwa-mem3/issues/131)) ([4d845c9](https://github.com/fg-labs/bwa-mem3/commit/4d845c92154d29f944191ec080009d916ce4b72b))

## [0.2.2](https://github.com/fg-labs/bwa-mem3/compare/v0.2.1...v0.2.2) (2026-06-08)


### Bug Fixes

* **lto:** pass explicit -flto=N on GCC to bypass jobserver under sandboxes ([#122](https://github.com/fg-labs/bwa-mem3/issues/122)) ([c6240a7](https://github.com/fg-labs/bwa-mem3/commit/c6240a74959e35a8fcbaa29a6f3b78c6a67ec7c3))
* **smem:** free lockstep SMEM caches at thread exit (closes [#116](https://github.com/fg-labs/bwa-mem3/issues/116)) ([#117](https://github.com/fg-labs/bwa-mem3/issues/117)) ([9454f10](https://github.com/fg-labs/bwa-mem3/commit/9454f106cba06ad3b139a01e5fc2b73a24356ce7))


### Performance

* **sort:** stabilize alnreg tie-breaks + drop in pdqsort at dedup-patch sort sites ([#123](https://github.com/fg-labs/bwa-mem3/issues/123)) ([85f8542](https://github.com/fg-labs/bwa-mem3/commit/85f8542080863c2fda9c198cb8743d18803e83e8))


### Documentation

* **bench:** inject generated divergence catalog + per-release concordance table ([#126](https://github.com/fg-labs/bwa-mem3/issues/126)) ([fea1c94](https://github.com/fg-labs/bwa-mem3/commit/fea1c9407a71b4e6c0efbd12aaa6a3a98b85861a))
* correct concordance claims and document supplementary-alignment divergence ([#125](https://github.com/fg-labs/bwa-mem3/issues/125)) ([8b2dc69](https://github.com/fg-labs/bwa-mem3/commit/8b2dc697c926c779b03be6369e8a08e95283b284))
* document bwa-mem3&lt;-&gt;bwa-mem2 non-bit-identity + auditable PR list ([#124](https://github.com/fg-labs/bwa-mem3/issues/124)) ([bffae5a](https://github.com/fg-labs/bwa-mem3/commit/bffae5a09267877fe514c458d4956b717bcefb8f))

## [0.2.1](https://github.com/fg-labs/bwa-mem3/compare/v0.2.0...v0.2.1) (2026-05-17)


### Bug Fixes

* **changelog:** strip preamble so release-please owns the file ([#112](https://github.com/fg-labs/bwa-mem3/issues/112)) ([56e580c](https://github.com/fg-labs/bwa-mem3/commit/56e580cf4f2515280440b2882a32c3a6a7b6d15c))
* **mapq:** propagate SMEM SA-count to seed n_hits so --supp-rep-hard-cap works ([#101](https://github.com/fg-labs/bwa-mem3/issues/101)) ([cca9d4f](https://github.com/fg-labs/bwa-mem3/commit/cca9d4f41023e501636d2b1ebd2d1e825f95e8e3))
* **smem:** track enc_qdb byte capacity separately from wsize_mem ([#100](https://github.com/fg-labs/bwa-mem3/issues/100)) ([ab922b6](https://github.com/fg-labs/bwa-mem3/commit/ab922b6af531ed4521b6cd9bd3c43a661d51a421))


### Documentation

* **readme:** add bioconda badges and install instructions ([#106](https://github.com/fg-labs/bwa-mem3/issues/106)) ([830276c](https://github.com/fg-labs/bwa-mem3/commit/830276ce01774805c20bfff69c63cfebab239166))
