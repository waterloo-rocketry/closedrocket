#include "_coder_GNC_codegen_SIL_info.h"
#include "emlrt.h"
#include "tmwtypes.h"

static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void);

static const mxArray *c_emlrtMexFcnResolvedFunctionsI(void) {
  const mxArray *nameCaptureInfo;
  const char_T *data[27] = {
      "789ced5d4d8ceb5615f6a0b67ae5f15e1fa5a50881a02c2ab1e8532699642645159afc4c"
      "32c9e4ff77f2f29af1d84ee2193b4e6c273399d5db8240884d570816"
      "4808a9124554482d749105a280ba60c182053fd2136a114845080958928ce349ec37263f"
      "b66fe2fbced9643c27f77ef71ee77e3939f7dc6362633fb14110c46d",
      "4291dfbfbc71f97a6b7c7d67fcfa31422b7afdc6f8d5a7bb56e549e2094d3b55ffadf12b"
      "25b464e65c562e5a24cf5cb5a4059e6d912d39df6f3384c84802d763"
      "e84b4d9de5983ccb33b9e98be4e88adf9b525d5d8c54a3bf834d863acd7579426c4a9311"
      "72d31757f6383298ef1333eca117bd3df4ef53f1ce97c453fbffcc0c",
      "3c553fb2b728701c23d62881661a4cabc60cffd3373b6f15f729c371281a4916bb943c99"
      "f7c0245ec5104fabbf17be1f7ca55a901851aaca222bc9d5904075f9"
      "e1dca56a8495a3dde3aac4f25d8e6d9dbe4c912d52a4a56aa34555c7f6aa1ad9ed2eaff4"
      "3fcb7ecfcc391ffdebe4fd372e5f5ff9fc1f3750e279a3ef7e15259e",
      "2aabc233bb0e3f6d807747a7cf52defdd079dc95eb9fb3b964b053da4e52d9c8641ce919"
      "38b3c641185ca3eaffc8a0bd5dfc893bde0f97c453fb8fcec053f5f7"
      "0acbf1243f2445aed626459297eef2a48c6afd7ee925b47c28bef54e13259e2ab8f3e1f9"
      "99bf9139f353e9e376b4123949bbfcbbfb4c181f3e1c18b477aa9f63",
      "f673f1dc8cf9a8faa9f13092cc0e994510a7c77164721c8bfaab6f9bc42b1ae269f516de"
      "c72bbb8dee212a1eb9f705b4bcfcd17b5f7909259e2ad8f3f2569e3b"
      "27b35b8d7299761fc738ff21c7378197d79697db26e7735b77ad9f8faa9f1a0f479ea1e3"
      "635ae81e73cc04efc726f1d2ba6b42f73e556fe1fd1bdaeb2e8f8a37",
      "9e46ccc3f1efbcbf8f124f15dc79b81327ebcde4d64976dfe5f3c4ebfe6c25158805f1e1"
      "e12383f64e8fb7b6c81edb20655668e1116fbd6f88a7d52fcb97137b"
      "558d4c77f995876abdbf8838bef0af177ef2359478aae0ce9f52b7ebbba0c242a2271533"
      "a54aa2d4167c31f0631db7aee79dcf0dddf5643e8aa6c750b59620ce",
      "f4c3ecf25b0726f1ecfedd3175bf788694ba22a3b451ed867adfeb7b9f43cbc35f7ffecf"
      "3750e2a9823b0f8783c1ba3fee8d0a312ed64ad70f05d9554a848087"
      "d79587cdc6133e3e633eaabe4dd2b53acbc98ca85cafca3f361b4f4819e269f516dcbf89"
      "c986770c577ff8072ffc14fc61c27a1ea6b3bd9d42b09f88d6b7d35b",
      "fd40863a49ef89517c78d8e9ebd82cefde9c317e553f02e784b336292959704edd57cb1b"
      "e269f5d6dcafea94d9144f18155fdcffdd1f90f2efc39b1f7851e2a9"
      "823bff9685a6ffa0cfb285f26e920ceccb814cd613db03fec5857f8dec7547f73a02675b"
      "d4707c02cf0c87e0d478c4a1219e566f21ff4e9b4d0d47601b8fb8ff",
      "0d88478c64decfe3f306787774fa54ce976ef47b72b25d0ad26eb1934b64188e001e5e17"
      "1e469567d6e99243e0d62818220a324f0ee77009bd2a3efe95493cda"
      "104fabb7e0be4d4c57bdce8a68fde357ff86d63fbe9979f1399478aa389597e7f58f430c"
      "d3f1c5ce36039cf7b059dc0a56e42ec96294efe0745e36eb1fdfd25d",
      "ebc7afea4956a44999ac91322f8c2214abca17192c89a7f68f70bf6e6cb2aac674a8f7eb"
      "7eb981d63f7efddf1fbe81124f15dc79982b45e55d575e8814723457"
      "0f31293e93cb629437013cfcffc76f170faf2a4e013c6c2f1ef0b022c0c38bf53f3068ff"
      "b8e64dcccbcb758e6d34e549ea8453f7efb286785abd05f74f6332d4",
      "7c8cfabcf2277ef37338af4cd8904f1cc8e5a2f983642195dbf1a73af57da975d1c168ff"
      "0e97f56cf6f3f0c919f350f5cc69bdd616196510b5e18809ed388e4c"
      "8e03b7ba3a53f76f68baea23e6c33dbf987fef4fb09f4758cfcb14bbeb0f309593646ffb"
      "3cd9dd4bb90edc111aa3b8f1c0a0bd53d7b5593f79defc629e949b35",
      "99e568468b7f64127fd1f885d9efd5b2219e566fc1fda3fb2d926729a93ab11dba73cba8"
      "f3dcfef38bbf409e1b617d7e45dcd5dfde2cbb859d8b1d4f2a1cca6e"
      "459ae56d021f3ec6c54f46c5c3a32f846e9b266587f370ce104fabb7f07b5431dbf8c6e1"
      "7adea3f7fa9b70de83b0de2f2ef4e221a19ecef8857dc11fe44a21fe",
      "7cb380915f8cfb7ab6eadcb3ead4399577edae9f769dffabfe81b27e1aea38f1b7377e06"
      "7162c27ade2d9de5b294cbe3f39c1e66a59ec7dd38487a248cce3be3"
      "b29e51c589a73264a7dd60c7e6199386785abd6d79c69a6f517cf38ce51f419ef148ace6"
      "e7cd20ddcfe699039ef2ba7cb94abbe28fd4798cf22a705ddf66e315",
      "cfe8aef5f352f553231995bd716cde5bcd104fabb7ed3e8e8b06e1cacffffd32f0f348ac"
      "e6e78b5497db77178abdbd3df745d1e78a5712c20e467916b8f2b3d9"
      "cfc7a766cc4bd54f8d64387a996d737dcd388e4c8e03b53f4d19e269f5b6dd4fd58a88e3"
      "cca8f7fb9eddf810f6fb08ebf7fb02d168f8fc3410e377e31d171523",
      "3d3b8d7e91c087af71897798f59fe7ad2f74993b200ab21ae970aaff8cb0be8536ef42b1"
      "1df6f18db7bff945f09f09ebfd675728c35562ac70de6b66b387645c"
      "ea7477b630daf71b18b477eaba36cbcb900fa77d1fe4c32d8707f9708a403edc62fd837f"
      "acfdbfd13c54bd8a5d231951209c7b5ee49e219e566fc77dbbb41de2",
      "7805eafc8cb752909f3112cbcf55ef77bddce969cb150d90914e5a124a54308f517ec6c0"
      "a0bd53d7352a5e867a1750ef621e3ca877a108d4bb58ac7f5cfc64a8"
      "cfb91c1ed4e7448307f5391581fa9c8bf58fcb79125479cd573efb094909c72cd91affdf"
      "a97e33c27cb9477ff7a83644982f873a9ef1aa0fea128dc4f27c6639",
      "51dc3c2d76fc7b7b671e5f36de0ed3c970001f5e1e18b477eaba7e60723e9f9d311f55ff"
      "48c69e8ea857c5d3bf3589776288a7d5db9cff387d6b2f0557defe7e"
      "f05de06dc2867a72bec3c304df144a2722b395e4dd8d10d3744780b7d795b7913da7713a"
      "163ec5d9abda277450dce3fafd041d53e3cad3df7d03f60b476279dc",
      "c35bdfd9aed0a218d92dee773abe482316cf6214f7c0757db74dce0bf60db5ef837d4373"
      "78b06fa808ec1b2ed6ffc0a0bd53fd68b3bc0ce750b4ef837328e6f0"
      "e01c8a22700e65b1fe0706ed1f575e867328daf7c13994e5f0e01c8a22700e65b1fe0706"
      "ed9dcac70f4cce67897d42c5e3837d42cbf2ed56b94f08f590ecc543",
      "e54767bcc502e565c3a47bb7d8c9677c5d3e447520beb1b6bc8d2a3f7aea37ba86b19d9a"
      "1f7d6c88a7d5db14e7786c781ae21d8a585ef7397150289793c1269b"
      "3b2989edf859d1d3e0a06edddaaf6f887b2c8607710f7bf020eea108c43d16eb7f60d0de"
      "a9feb3593ebeadbbd6cf47d54ffd52675b3de7ee0fbe6688a7d5db16",
      "e7185a0fe57352a06ea8bd78a8fc655fb6dbc98baed00e4b7bce3242b0eee1b752508763"
      "6d7919ce7d2f8707e7bed1e0c1b96f45e0dcf762fde372eedbacdf0c"
      "cf47d1be0f9e8f620e0ff6031581e7a32cd63f3c5f5b91a767cc43d5f30c29d58ec9cb5a"
      "a2d3f8474be2ebc5085f15abf81861ddc191c9ba22a3b4b9b21fe27a",
      "a2a8cf9f3cf5eb8fe0fc0961031fb319dfa94b94dc07799a89e56391ad0c2d47818fd78d"
      "8fcd7e1e9e9d310f557fc52716e7d1e9655e5e7650fce27a5e5ed1f9"
      "6dd8ffb3170f153fa74936594ad00d0f158fb543debca72976388cf2e770e167b3fef2bc"
      "e7024701144a10456654e4c242fc79e327f56e8b1a59a3d6245b34c7",
      "58761ff333f055bd8571a8b119d1f2f2ed27d1facd7f7fed9f1d9478aa389597e7cdcb10"
      "63f5134fa8e9e22f28b27ed8edf502b9fc2601bcbc6eebd9a9718c55"
      "c595218e612f1ec431148138c662fd0f0cda3b755d9be5e55bba6bfd7c543dd43782fa46"
      "f3e0417d2345a0bed162fde3e2273b3dbebc2c2f437cf9fa79417c59",
      "11dcf919f7f832aeeb1bd5e743e304aea09eb3d57c4d19e269f576f9d1fa839cb8f2b5b0"
      "f710f89ab07e3d0632e7c2b15c38cbc9f9cd462f17cc6f735205f299"
      "d7ce9f361bdf807c66edfb209fd91c1ee4332b02f9cc8bf58f4b7ec691c979dcd05d4fe6"
      "a1682efd759e6c38357f19617cf9d1df3943bbe11e5f7ef3ad7f407c",
      "99b0e17c9fdb7db029bbb36d36722ad25189f61eecc7c01f5e3bfe45f55c57954fac7eae"
      "ab5ee6e56507c52baee5e5c7255ef1f0e60710af20ace7e756982fe7"
      "64399d4efbf2ee42afdcebc68b71f08fd72e5e7164721ea8fce355c527c03fb6170ffc63"
      "45c03f5eacff81417ba7ae6ba853b41c1ed429428307758a14813a45",
      "8bf58f8b9fecf438c6aaf819e21868f0208ea108c43116eb1fd7f5fdc0e4bce039544678"
      "5a3d3c87ca1a3cc8d350049e43b558ff0383f6f3daf1be41ff77747a"
      "0bd6f9e4cfa18b4d330da6551b3617fb974b1bce072e8707e703edc583f3818ac0f9c0c5"
      "fa1f18b47f5c7979debc66955cc635fd91f1b2248bdd516925abee5f",
      "d5104fabb78197c7a6431b8f46cdcbe5f7ff0abc4c58cfcbdeeec57186ec35c2d17cee24"
      "c5a6530c9d3b0c002faf2b2fa33abf3df5cb5c16863fce394ef9bf53"
      "e3d2089feb787d7c636c45dce31a89fb10d71889d53c9ded784b89403dec4ff6fdfd9d70"
      "9ccfb9f91c063cfd3f82f44909",
      ""};
  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(&data[0], 77296U, &nameCaptureInfo);
  return nameCaptureInfo;
}

mxArray *emlrtMexFcnProperties(void) {
  mxArray *xEntryPoints;
  mxArray *xInputs;
  mxArray *xResult;
  const char_T *epFieldName[7] = {
      "QualifiedName",    "NumberOfInputs", "NumberOfOutputs", "ConstantInputs",
      "ResolvedFilePath", "TimeStamp",      "Visible"};
  const char_T *propFieldName[7] = {
      "Version",      "ResolvedFunctions", "Checksum", "EntryPoints",
      "CoverageInfo", "IsPolymorphic",     "AuxData"};
  uint8_T v[216] = {
      0U,   1U,   73U,  77U,  0U,   0U,   0U,   0U,   14U,  0U,   0U,   0U,
      200U, 0U,   0U,   0U,   6U,   0U,   0U,   0U,   8U,   0U,   0U,   0U,
      2U,   0U,   0U,   0U,   0U,   0U,   0U,   0U,   5U,   0U,   0U,   0U,
      8U,   0U,   0U,   0U,   1U,   0U,   0U,   0U,   1U,   0U,   0U,   0U,
      1U,   0U,   0U,   0U,   0U,   0U,   0U,   0U,   5U,   0U,   4U,   0U,
      17U,  0U,   0U,   0U,   1U,   0U,   0U,   0U,   17U,  0U,   0U,   0U,
      67U,  108U, 97U,  115U, 115U, 69U,  110U, 116U, 114U, 121U, 80U,  111U,
      105U, 110U, 116U, 115U, 0U,   0U,   0U,   0U,   0U,   0U,   0U,   0U,
      14U,  0U,   0U,   0U,   112U, 0U,   0U,   0U,   6U,   0U,   0U,   0U,
      8U,   0U,   0U,   0U,   2U,   0U,   0U,   0U,   0U,   0U,   0U,   0U,
      5U,   0U,   0U,   0U,   8U,   0U,   0U,   0U,   1U,   0U,   0U,   0U,
      0U,   0U,   0U,   0U,   1U,   0U,   0U,   0U,   0U,   0U,   0U,   0U,
      5U,   0U,   4U,   0U,   14U,  0U,   0U,   0U,   1U,   0U,   0U,   0U,
      56U,  0U,   0U,   0U,   81U,  117U, 97U,  108U, 105U, 102U, 105U, 101U,
      100U, 78U,  97U,  109U, 101U, 0U,   77U,  101U, 116U, 104U, 111U, 100U,
      115U, 0U,   0U,   0U,   0U,   0U,   0U,   0U,   80U,  114U, 111U, 112U,
      101U, 114U, 116U, 105U, 101U, 115U, 0U,   0U,   0U,   0U,   72U,  97U,
      110U, 100U, 108U, 101U, 0U,   0U,   0U,   0U,   0U,   0U,   0U,   0U};
  xEntryPoints =
      emlrtCreateStructMatrix(1, 2, 7, (const char_T **)&epFieldName[0]);
  xInputs = emlrtCreateLogicalMatrix(1, 6);
  emlrtSetField(xEntryPoints, 0, "QualifiedName",
                emlrtMxCreateString("controller_codegen_entry"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(6.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(4.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(
      xEntryPoints, 0, "ResolvedFilePath",
      emlrtMxCreateString("C:\\Users\\trist\\Documents\\GitHub\\simulink-"
                          "canards\\gnc\\control\\controller_codegen_entry.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp",
                emlrtMxCreateDoubleScalar(740189.6711805556));
  emlrtSetField(xEntryPoints, 0, "Visible", emlrtMxCreateLogicalScalar(true));
  xInputs = emlrtCreateLogicalMatrix(1, 7);
  emlrtSetField(xEntryPoints, 1, "QualifiedName",
                emlrtMxCreateString("navigation_codegen_entry"));
  emlrtSetField(xEntryPoints, 1, "NumberOfInputs",
                emlrtMxCreateDoubleScalar(7.0));
  emlrtSetField(xEntryPoints, 1, "NumberOfOutputs",
                emlrtMxCreateDoubleScalar(8.0));
  emlrtSetField(xEntryPoints, 1, "ConstantInputs", xInputs);
  emlrtSetField(xEntryPoints, 1, "ResolvedFilePath",
                emlrtMxCreateString(
                    "C:\\Users\\trist\\Documents\\GitHub\\simulink-"
                    "canards\\gnc\\navigation\\navigation_codegen_entry.m"));
  emlrtSetField(xEntryPoints, 1, "TimeStamp",
                emlrtMxCreateDoubleScalar(740201.00158564816));
  emlrtSetField(xEntryPoints, 1, "Visible", emlrtMxCreateLogicalScalar(true));
  xResult =
      emlrtCreateStructMatrix(1, 1, 7, (const char_T **)&propFieldName[0]);
  emlrtSetField(xResult, 0, "Version",
                emlrtMxCreateString("25.2.0.2998904 (R2025b)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions",
                (mxArray *)c_emlrtMexFcnResolvedFunctionsI());
  emlrtSetField(xResult, 0, "Checksum",
                emlrtMxCreateString("bKrbaiBFSs4CbFe2poCkxG"));
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  emlrtSetField(xResult, 0, "AuxData",
                emlrtMxCreateRowVectorUINT8((const uint8_T *)&v, 216U));
  return xResult;
}
